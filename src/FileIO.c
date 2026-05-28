#include "FileIO.h"
#include "Atoms.h"
#include "ErrorM.h"
#include "Input.h"
#include "InputXYZ.h"
#include "Utils.h"
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const size_t BUFFER_SIZE = 256;
const size_t ARR_BUFFER_SIZE = 32;

static void fopen_error(char *filename, FILE *file, char *base_msg)
{
    if (file == NULL) {
        fprintf(stderr, "%s%s: %s\n", base_msg, filename, strerror(errno));
        clean_and_error(errno);
    }
    return;
}

/**
 * @brief safely logs a formatted message to the given stream, ensuring that the entire message is
 * written; exits with an error if formatting or writing fails
 *
 * @param stream file stream to write to
 * @param fmt printf-like format string
 * @param ... arguments for format string
 */
void safe_log(FILE *stream, const char *fmt, ...)
{
    // write formatted message to buffer first to ensure atomicity of log lines and minimize chance
    // of partial writes
    char buffer[1024];
    va_list args;
    va_start(args, fmt);
    int n = vsnprintf(buffer, sizeof(buffer), fmt, args);
    va_end(args);

    // check for formatting errors before writing
    if (n < 0) {
        fprintf(stderr, "Formatting error in safe_log, format %s: %s\n", fmt, strerror(errno));
        clean_and_error(errno);
    }
    if (n >= (int)sizeof(buffer)) {
        fprintf(stderr, "Log line too long (%d chars), buffer size %zu\n", n, sizeof(buffer));
        clean_and_error(EXIT_FAILURE);
    }

    // write buffer to stream and check for write errors
    size_t written = fwrite(buffer, 1, (size_t)n, stream);
    if (written != (size_t)n) {
        // wrote only part of the line
        // don't error for now
        fprintf(stderr,
                "Failed to write complete log line to stream - wrote %zu of %d characters: %s\n",
                written, n, strerror(errno));
        clean_and_error(errno);
    }

    int fres = fflush(stream);
    if (fres == EOF) {
        perror("Failed to flush log");
        clean_and_error(errno);
    }
}

void open_log_files(struct LoggingState *ls, unsigned flavor)
{
    for (int i = 0; i < ls->out_formats_cnt; i) {
        OutputFormat *format = &(ls->out_formats[i]);
        if (format->type == OUTPUT_FORMAT_CSV) {
            format->csv.file = fopen(format->csv.filename, "w+");
            fopen_error(format->csv.filename, format->csv.file, "Failed to open csv file, ");
            output_csv_header(format->csv.file, format);
        } else if (format->type == OUTPUT_FORMAT_STEPS_CSV) {
            format->csv.file = fopen(format->steps.filename, "w+");
            fopen_error(format->steps.filename, format->csv.file,
                        "Failed to open steps csv file, ");

            if (flavor == FLAVOR_KMC) {
                output_kmc_steps_header(format->steps.file, format->steps.coordination);
            } else if (flavor == FLAVOR_MC) {
                output_mc_steps_header(format->steps.file, format->steps.coordination);
            }
        }
    }
}

void simulation_parameters_from_file(char *filename, struct SimulationConfig *inputs,
                                     struct LoggingState *ls)
{
    FILE *input_file = fopen(filename, "r");
    fopen_error(filename, input_file, "Failed to open input file, ");

    process_in_file(input_file, inputs, ls);
    fclose(input_file);
    return;
}

void process_in_file(FILE *input_file, struct SimulationConfig *inputs, struct LoggingState *ls)
{
    ParseContext ctx = {0};
    parse_input_file(input_file, &ctx, inputs, ls);
    finalize_config(&ctx, inputs, ls);
    clean_ctx(&ctx);
    return;
}

/**
 * @brief parses an .xyz file to populate the simulation state and environment
 *
 * @param input_file
 * @param ss
 * @param se
 * @param ls
 * @return true successfully parsed xyz file and populated simulation state and env
 * @return false failed to parse xyz file
 */
bool process_xyz_file(FILE *input_file, struct SimulationState *ss, struct SimulationEnv *se,
                      struct LoggingState *ls)
{
    // processes file with .xyz format (number of atoms / comment / type x y z)
    // LoggingState only needed for xyz_framenum

    char command_string[BUFFER_SIZE];
    char *ptr; // for fgets return values

    ptr = fgets(command_string, (int)BUFFER_SIZE, input_file);

    if (!ptr) {
        fprintf(stderr, "Empty file or read error\n");
        return false;
    }

    // first line is the number of atoms (lines with atom info)
    char *endptr = NULL;
    long natoms_long = strtol(command_string, &endptr, 10);
    if (endptr == command_string || natoms_long <= 0) {
        fprintf(stderr, "First line does not contain a valid atom count\n");
        return false;
    }
    int nremain = (int)natoms_long;

    // read comment line
    unsigned int comment_buffer_multiplier = 3;
    char comment_string[comment_buffer_multiplier * BUFFER_SIZE];
    ptr = fgets(comment_string, (int)(comment_buffer_multiplier * BUFFER_SIZE), input_file);

    if (!ptr) {
        fprintf(stderr, "Error on reading the comment/header line\n");
        return false;
    }

    // TODO: provide flag to control whether comments get parsed
    // act upon key-value pairs - simulation variables, Properties
    // capture any simulation variables in the comment line that are used to continue a started
    // simulation ss->elapsed_stime, ss->temperature, ss->overpotential; all doubles ls->framenum;
    // int ss->iter; unsigned long
    size_t kvpairs_cnt = 0;
    struct KV *kvpairs = NULL; // pointer to array of KVs
    int pc = parse_comment(comment_string, &kvpairs, &kvpairs_cnt);
    if (pc) {
        return pc;
    }

    int has_props = 0;
    struct KV kv;
    PropertyDesc *properties = NULL;
    int properties_cnt;
    for (int i = 0; i < (int)kvpairs_cnt; i++) {
        kv = kvpairs[i];
        if (strncmp(kv.key, "properties", 10) == 0) {
            int pp = parse_properties_value(kv.value, &properties, &properties_cnt);
            if (pp) {
                clean_xyz_structs(kvpairs, kvpairs_cnt, properties);
                return false;
            }
            has_props = 1;
        } else if (strncmp(kv.key, "time", 4) == 0) {
            ss->elapsed_stime = strtod(kv.value, NULL);
        } else if (strncmp(kv.key, "temperature", 11) == 0) {
            ss->temperature = strtod(kv.value, NULL);
        } else if (strncmp(kv.key, "potential", 9) == 0) {
            ss->overpotential = strtod(kv.value, NULL);
        } else if (strncmp(kv.key, "iteration", 9) == 0) {
            ss->iter = strtoul(kv.value, NULL, 10);
        } else if (strncmp(kv.key, "energy", 6) == 0) {
            ss->total_internal_energy = strtod(kv.value, NULL);
        }
    }

    char *tokens[64]; // array of char arrays (array of char pointers)
    for (int i = 0; i < nremain; i++) {
        ptr = fgets(command_string, (int)BUFFER_SIZE, input_file);

        // if EOF or failure
        if (!ptr) {
            fprintf(stderr, "Input parsing failed - Ran into EOF, expected %d atoms remaining\n",
                    nremain);
            clean_xyz_structs(kvpairs, kvpairs_cnt, properties);
            return false;
        }

        int ntok = tokenize_line(command_string, tokens,
                                 (int)(sizeof(tokens) / sizeof(tokens[0]))); // aka 64

        Atom temp_atom = {0};
        if (has_props) {
            // Check we have enough tokens for declared properties
            int ft = fill_atom_from_tokens(&temp_atom, tokens, ntok, se->atom_names,
                                           se->atom_names_cnt, properties, properties_cnt);
            if (ft) {
                fprintf(stderr, "Atom line %d has %d tokens, expected >= %d\n", i, ntok, ft);
            }
        } else {
            fill_atom_from_xyz(&temp_atom, tokens, ntok, se->atom_names, se->atom_names_cnt);
        }
        cartesian2lattice_site(temp_atom.cartesian, se->invert_primitive_basis, temp_atom.lattice);
        add_atom(temp_atom.lattice[0], temp_atom.lattice[1], temp_atom.lattice[2], temp_atom.type,
                 NORMAL, ss, se);
    }

    ptr = fgets(command_string, (int)BUFFER_SIZE, input_file);
    // end of file is not reached
    // or a final blank line is not last line
    if (!feof(input_file)) {
        int nonspace = 0;
        int clean = 0;
        for (size_t i = 0; i < strlen(command_string); i++) {
            if (!isspace(command_string[i])) {
                nonspace++;
            }
        }
        if (nonspace) {
            clean = 1;
        } else {
            ptr = fgets(command_string, (int)BUFFER_SIZE, input_file);
            if (!feof(input_file)) {
                clean = 1;
            }
        }
        if (clean) {
            fprintf(stderr,
                    "Unparsed content remains in input file after %d reads; first line was not "
                    "accuraten\n",
                    nremain);
            clean_xyz_structs(kvpairs, kvpairs_cnt, properties);
            return false;
        }
    }

    clean_xyz_structs(kvpairs, kvpairs_cnt, properties);

    safe_log(stdout, "Successfully read %ld atoms from .xyz file\n", ss->atom_cnt);
    return true;
}

/* === CSV fields mapping === */

const char *fstring_to_buffer(const char *fmt, ...)
{
    va_list vargs;
    va_start(vargs, fmt);
    int len = vsnprintf(NULL, 0, fmt, vargs);
    va_end(vargs);

    // len is negative when vsnprintf errors
    // but errors may depend on implementation
    if (len < 0) {
        return NULL;
    }

    char *buf = (char *)malloc((size_t)len + 1);
    if (!buf) {
        return NULL;
    }

    va_start(vargs, fmt);
    vsnprintf(buf, (size_t)len + 1, fmt, vargs);
    va_end(vargs);

    return buf;
}

static const char *get_iteration(const struct SimulationState *ss, const struct CsvLsView *view)
{
    (void)view;
    return fstring_to_buffer("%lu", ss->iter);
}

const char *get_time(const struct SimulationState *ss, const struct CsvLsView *view)
{
    return fstring_to_buffer("%.*le", view->precision, ss->elapsed_stime);
}

const char *get_energy(const struct SimulationState *ss, const struct CsvLsView *view)
{
    (void)view;
    return fstring_to_buffer("%lf", ss->total_internal_energy);
}

const char *get_temperature(const struct SimulationState *ss, const struct CsvLsView *view)
{
    (void)view;
    return fstring_to_buffer("%lf", ss->temperature);
}

const char *get_overpotential(const struct SimulationState *ss, const struct CsvLsView *view)
{
    return fstring_to_buffer("%.*le", view->precision, ss->overpotential);
}

const char *get_atoms(const struct SimulationState *ss, const struct CsvLsView *view)
{
    (void)view;
    return fstring_to_buffer("%ld", ss->atom_cnt);
}

const char *get_atoms_dissolved(const struct SimulationState *ss, const struct CsvLsView *view)
{
    (void)view;
    return fstring_to_buffer("%d", ss->total_atoms_dissolved);
}

const char *get_mmc_steps(const struct SimulationState *ss, const struct CsvLsView *view)
{
    (void)view;
    return fstring_to_buffer("%lu", ss->mmc_steps);
}

const CsvFieldFunc csv_field_map[] = {{"iter", get_iteration, FLAVOR_UNDEFINED},
                                      {"time", get_time, FLAVOR_UNDEFINED},
                                      {"energy", get_energy, FLAVOR_UNDEFINED},
                                      {"temperature", get_temperature, FLAVOR_UNDEFINED},
                                      {"overpotential", get_overpotential, FLAVOR_KMC},
                                      {"atoms", get_atoms, FLAVOR_UNDEFINED},
                                      {"atoms_dissolved", get_atoms_dissolved, FLAVOR_KMC},
                                      {"mmc_steps", get_mmc_steps, FLAVOR_MC}};

const size_t CSV_FIELD_FUNCS_COUNT = sizeof(csv_field_map) / sizeof(CsvFieldFunc);

static char *schedule_list_to_string(OutputSchedule *schedule)
{
    char *buf = NULL;
    // max size of each item in list
    size_t item_size = 32;
    if (schedule->list != NULL) {
        buf = (char *)malloc((size_t)schedule->list_len * item_size * sizeof(char));
        buf[0] = '\0';
        int printed_size = snprintf(buf, item_size, "%lf ", schedule->list[0]);
        // printed_size does not include null terminating character
        if (printed_size >= (int)item_size) {
            fprintf(stderr, "Error - Output schedule list item string too long (>%zu chars): %lf\n",
                    item_size - 2, schedule->list[0]);
            clean_and_error(EXIT_FAILURE);
        }
        for (int i = 1; i < schedule->list_len; i++) {
            char item[item_size];
            strcat(buf, ", ");
            printed_size = snprintf(item, item_size, "%lf ", schedule->list[i]);
            if (printed_size >= (int)item_size) {
                fprintf(stderr,
                        "Error - Output schedule list item string too long (>%zu chars): %lf\n",
                        item_size - 2, schedule->list[0]);
                clean_and_error(EXIT_FAILURE);
            }
            strcat(buf, item);
        }
    }
    return buf;
}

// print a lot of information to the log
// use simulation structs when possible to show values have been set, otherwise grab from inputs
void input_logging(struct SimulationConfig *inputs, struct SimulationState *ss,
                   struct SimulationEnv *se, struct LoggingState *ls)
{
    safe_log(ls->sim_log, "Successfully read input file and preprocessed\n");
    safe_log(ls->sim_log, "System size is %d x %d x %d\n", se->system_size_x, se->system_size_y,
             se->system_size_z);

    safe_log(ls->sim_log, "Crystal structure is ");
    switch (inputs->lattice_type) {
    case FCC:
        safe_log(ls->sim_log, "FCC");
        break;
    case BCC:
        safe_log(ls->sim_log, "BCC");
        break;
    case SC:
        safe_log(ls->sim_log, "SC");
        break;
    }
    safe_log(ls->sim_log, "\n");

    safe_log(ls->sim_log, "Initializing atom types: ");
    for (int i = 0; i < se->num_elements; i++) {
        safe_log(ls->sim_log, "%s ", se->atom_names[i]);
    }
    safe_log(ls->sim_log, "\nComposition: ");
    for (int i = 0; i < se->num_elements; i++) {
        safe_log(ls->sim_log, "%lf ", se->substrate_composition[i]);
    }

    safe_log(ls->sim_log, "\nSolubility: ");
    for (int i = 0; i < se->num_elements; i++) {
        safe_log(ls->sim_log, "%s ", se->is_soluble[i] ? "true" : "false");
    }

    safe_log(ls->sim_log, "\nBond energies\n");
    int bond_idx, env_idx;
    for (int nn_level = 0; nn_level < se->num_nn_levels; nn_level++) {
        for (int elem_a = 0; elem_a < se->num_elements; elem_a++) {
            for (int elem_b = elem_a; elem_b < se->num_elements; elem_b++) {
                bond_idx = get_bond_index(elem_a, elem_b, se->num_elements);
                env_idx = nn_bondidx_2_envidx(nn_level, bond_idx, se->num_bond_types);
                safe_log(ls->sim_log, "%s-%s: %lf\n", se->atom_names[elem_a],
                         se->atom_names[elem_b], se->nn_energy[env_idx]);
            }
        }
    }
    char flavor_str[4];
    switch (inputs->flavor) {
    case FLAVOR_KMC:
        strncpy(flavor_str, "KMC", sizeof(flavor_str));
        break;
    case FLAVOR_MC:
        strncpy(flavor_str, "MC", sizeof(flavor_str));
        break;
    default:
        strncpy(flavor_str, "IDK", sizeof(flavor_str));
        break;
    }
    safe_log(ls->sim_log, "Flavor is %s\n", flavor_str);

    safe_log(ls->sim_log, "Temperature is %lf K\n", ss->temperature);

    if (se->overpotential_ramp_rate > 0.) {
        safe_log(ls->sim_log, "Potential sweep from %lfeV to %lfeV at %lfeV/s\n", ss->overpotential,
                 se->overpotential_ramp_rate, se->max_overpotential);
    } else {
        safe_log(ls->sim_log, "Potential constant [eV] at %lf\n", ss->overpotential);
    }

    for (int i = 0; i < ls->out_formats_cnt; i++) {
        OutputFormat *format = &(ls->out_formats[i]);
        switch (format->type) {
        case OUTPUT_FORMAT_CSV:
            safe_log(ls->sim_log, "CSV output to file %s\n", format->csv.filename);
            break;
        case OUTPUT_FORMAT_STEPS_CSV:
            safe_log(ls->sim_log, "Steps CSV output to file %s\n", format->steps.filename);
            break;
        case OUTPUT_FORMAT_XYZ:
            safe_log(ls->sim_log, "XYZ output to file prefix %s\n", format->xyz.prefix);
            break;
        }

        if (format->type == OUTPUT_FORMAT_CSV || format->type == OUTPUT_FORMAT_XYZ) {
            OutputSchedule schedule =
                (format->type == OUTPUT_FORMAT_CSV) ? format->csv.schedule : format->xyz.schedule;
            char *buf = schedule_list_to_string(&schedule);

            switch (schedule.mode) {
            case OUTPUT_SCHEDULE_INTERVAL_ITERATION:
                safe_log(
                    ls->sim_log,
                    "Recording data at linear intervals [iterations] from at %.2lf increments\n",
                    schedule.interval);
                break;
            case OUTPUT_SCHEDULE_INTERVAL_TIME:
                safe_log(ls->sim_log,
                         "Recording data at linear intervals [s] from at %.2lg multiples\n",
                         schedule.interval);
                break;
            case OUTPUT_SCHEDULE_LIST_ITERATION:
                safe_log(ls->sim_log, "Recording data at iterations %s\n", buf);
                break;
            case OUTPUT_SCHEDULE_LIST_TIME:
                safe_log(ls->sim_log, "Recording data at times %s\n", buf);
                break;
            default:
                fprintf(stderr, "Invalid output schedule mode\n");
                clean_and_error(EXIT_FAILURE);
                break;
            }
        }
    }

    safe_log(ls->sim_log, "Random seed is %u\n", se->rand_seed);

    switch (inputs->geometry) {
    case GEOMETRY_FLAT_SHEET:
        safe_log(ls->sim_log, "Initialized flat sheet with monolayer depth %d\n",
                 inputs->geometry_param);
        break;
    case GEOMETRY_CLUSTER:
        safe_log(ls->sim_log, "Initialized spherical cluster with radius %d\n",
                 inputs->geometry_param);
        break;
    case GEOMETRY_FROM_FILE:
        safe_log(ls->sim_log, "Initialized user-defined structure with filename %s\n",
                 inputs->atoms_filename);
        break;
    }

    safe_log(ls->sim_log, "Atoms created, %ld total\n", ss->atom_cnt);
}

/* === Data output === */

/**
 * @brief Get the precision that should be used for printing @p total with @c e formatting to
 * resolve the @p increment plus increment precision @p incr_precision. Assumes the increment
 * has already been added to the total.
 *
 * @param total
 * @param increment
 * @param incr_precision
 * @return int Resolving precision. Returns -1 if an error occurs.
 */
int get_precision(double total, double increment, int incr_precision)
{
    // precision is the number of digits to appear after the decimal point
    // cast truncates value
    if ((total == 0) || (increment == 0)) {
        return incr_precision;
    }
    double ltotal = log10(fabs(total));
    double linc = log10(fabs(increment));
    if (isnan(ltotal) || isnan(linc)) {
        return -1;
    }
    int log_diff = (int)fmax((int)ltotal - (int)linc, 0);
    // subtract one for the digit in front of the decimal
    int precision = log_diff + incr_precision;
    return precision;
}

bool output_log_file(FILE *sim_log, int frame_num, unsigned long int iter, double elapsed_stime,
                     double temperature, double overpotential, long int atom_cnt,
                     double total_internal_energy)
{
    safe_log(sim_log, "![%d]\t", frame_num);
    safe_log(sim_log, "iteration = %lu\t", iter);
    safe_log(sim_log, "time = %le [s]\ttemperature = %lf [K]\tpotential = %lf [eV]\t",
             elapsed_stime, temperature, overpotential);
    safe_log(sim_log, "atoms = %ld\tinternal energy = %lf [eV]\n", atom_cnt, total_internal_energy);
    return true;
}

// ENHANCE: make fields of union own type, so it can passed into this
void output_csv_header(FILE *csv_file, OutputFormat *format)
{
    safe_log(csv_file, "frame,");
    for (int i = 0; i < format->csv.field_count; ++i) {
        safe_log(csv_file, "%s%s", format->csv.field_names[i],
                 (i < format->csv.field_count - 1) ? "," : "\n");
    }
}

void log_state_csv(OutputFormat *format, double stime_precision, double overpot_precision,
                   struct SimulationState *ss)
{
    safe_log(format->csv.file, "%d,", format->csv.frame_num);
    for (int i = 0; i < format->csv.field_count; ++i) {
        // malloc'd string, needs to be free'd when done
        struct CsvLsView *view = malloc(sizeof(struct CsvLsView));
        if (strcmp(format->csv.field_names[i], "time") == 0) {
            view->precision = stime_precision;
        } else if (strcmp(format->csv.field_names[i], "overpotential") == 0) {
            view->precision = overpot_precision;
        }
        const char *value_str = format->csv.field_funcs[i](ss, view);
        if (!value_str) {
            fprintf(stderr, "Error formatting csv field %s: %s\n", format->csv.field_names[i],
                    strerror(errno));
            free((void *)value_str);
            clean_and_error(EXIT_FAILURE);
        }
        safe_log(format->csv.file, "%s%s", value_str,
                 (i < format->csv.field_count - 1) ? "," : "\n");
        free((void *)value_str);
    }
}

void output_kmc_steps_header(FILE *csv_file, const bool output_coord)
{
    safe_log(csv_file, "iter,sim_time,energy,u1,v1,w1,u2,v2,w2");
    if (output_coord) {
        safe_log(csv_file, ",coordination");
    }
    safe_log(csv_file, "\n");
}

// output to csv:
// iteration number, simulation time, system energy (per atom?), x1, y1, z1, x2, y2, z2
// and atom ids at some point
void log_kmc_steps(FILE *csv_file, const StepData *step_data, double sim_time_precision)
{
    safe_log(csv_file, "%lu,", step_data->iter);
    safe_log(csv_file, "%.*le,", sim_time_precision, step_data->kmc.sim_time);
    safe_log(csv_file, "%lf,", step_data->sys_energy);
    safe_log(csv_file, "%d,%d,%d,", step_data->uvw1[0], step_data->uvw1[1], step_data->uvw1[2]);
    if (!step_data->kmc.is_evap) {
        safe_log(csv_file, "%d,%d,%d", step_data->uvw2[0], step_data->uvw2[1], step_data->uvw2[2]);
    } else {
        safe_log(csv_file, ",,");
    }
    if (step_data->coord > 0) {
        safe_log(csv_file, ",%d", step_data->coord);
    }
    safe_log(csv_file, "\n");
}

void output_mc_steps_header(FILE *csv_file, bool output_coord)
{
    safe_log(csv_file, "iter,energy,deltaE,performed,u1,v1,w1,u2,v2,w2");
    if (output_coord) {
        safe_log(csv_file, ",coordination");
    }
    safe_log(csv_file, "\n");
}

// output to csv:
// MCSS, system energy (per atom?), uvw1, uvw2
void log_mc_steps(FILE *csv_file, const StepData *step_data)
{
    safe_log(csv_file, "%lu,", step_data->iter);
    safe_log(csv_file, "%lf,", step_data->sys_energy);
    safe_log(csv_file, "%lf,", step_data->mc.deltaE);
    safe_log(csv_file, "%d,", step_data->mc.performed);
    safe_log(csv_file, "%d,%d,%d,", step_data->uvw1[0], step_data->uvw1[1], step_data->uvw1[2]);
    safe_log(csv_file, "%d,%d,%d", step_data->uvw2[0], step_data->uvw2[1], step_data->uvw2[2]);
    if (step_data->coord > 0) {
        safe_log(csv_file, ",%d", step_data->coord);
    }
    safe_log(csv_file, "\n");
}

bool write_xyz_file(char *xyz_prefix, int frame_num, char *suffix, int stripped,
                    struct SimulationState *ss, struct SimulationEnv *se)
{
    bool is_extended = 1;

    char filename_full[BUFFER_SIZE];
    int n = snprintf(filename_full, BUFFER_SIZE, "%s_%d_%s.xyz", xyz_prefix, frame_num, suffix);
    if ((size_t)n >= BUFFER_SIZE) {
        fprintf(stderr, "Error - Output filename too long (>%zu): %s_%d_%s.xyz\n", BUFFER_SIZE,
                xyz_prefix, frame_num, suffix);
        clean_and_error(EXIT_FAILURE);
    }
    FILE *file = fopen(filename_full, "w+");
    if (file == NULL) {
        printf("Error - Couldn't open output file %s\n", filename_full);
        fprintf(stderr, "Couldn't open file %s: %s\n", filename_full, strerror(errno));
        clean_and_error(errno);
    }

    /* format:
            [number of atoms]
            [comment line - exactly one line]
            [element] [x] [y] [z]
    */

    // identify atoms to print before writing
    char *is_undercoord = (char *)malloc((size_t)ss->atom_cnt * sizeof(char));
    long strip_cnt = 0;

    if (stripped) {
        for (long i = 0; i < ss->atom_cnt; ++i) {
            int coord = get_coordination(i, ss, se);
            is_undercoord[i] = (char)(coord < se->num_transition_vectors);
            strip_cnt += is_undercoord[i];
        }
    }

    // start with number of atoms
    safe_log(file, "%ld\n", stripped ? strip_cnt : ss->atom_cnt);

    if (is_extended) {
        // using extended XYZ format
        // https://docs.ovito.org/reference/file_formats/input/xyz.html#file-formats-input-xyz-extended-format

        // 3x3 matrix - rows are cell vectors [preferred]
        safe_log(file, "Lattice=\"%lf %lf %lf %lf %lf %lf %lf %lf %lf\" ",
                 se->simbox_vectors_cart[0][0], se->simbox_vectors_cart[0][1],
                 se->simbox_vectors_cart[0][2], se->simbox_vectors_cart[1][0],
                 se->simbox_vectors_cart[1][1], se->simbox_vectors_cart[1][2],
                 se->simbox_vectors_cart[2][0], se->simbox_vectors_cart[2][1],
                 se->simbox_vectors_cart[2][2]);

        safe_log(file, "Origin=\"%lf %lf %lf\" ", se->simbox_origin_cart[0],
                 se->simbox_origin_cart[1], se->simbox_origin_cart[2]);
        safe_log(file, "pbc=\"T T T\" ");
        safe_log(file, "Properties=id:I:1:species:S:1:pos:R:3 ");
    }
    safe_log(file, "frame=%d iteration=%lu time=%le temperature=%lf potential=%lf energy=%lf\n",
             frame_num, ss->iter, ss->elapsed_stime, ss->temperature, ss->overpotential,
             ss->total_internal_energy);

    Atom **atoms = ss->atom_arr;
    for (long i = 0; i < ss->atom_cnt; ++i) {
        // if stripped, only output atoms with < max coordination
        // what is max coordination?
        int print = stripped ? is_undercoord[i] : 1;
        if (print) {
            safe_log(file, "%d %s %lf %lf %lf\n", i, se->atom_names[atoms[i]->type],
                     atoms[i]->cartesian[0], atoms[i]->cartesian[1], atoms[i]->cartesian[2]);
        }
    }

    fclose(file);

#if (!defined(NDEBUG)) && defined(DUMP_CORE)
    fprintf(stderr, "Creating core dump for frame %d\n", frame_num);
    fflush(stderr);
    create_coredump();
#endif

    return true;
}

static int check_and_advance_checkpoint(OutputSchedule *sched, bool *state,
                                        struct SimulationState *ss)
{
    OutputScheduleMode mode = sched->mode;
    int iter_reached = fabs(sched->next_checkpoint - (double)ss->iter) < FABS_TOL;
    int stime_reached = ss->elapsed_stime >= sched->next_checkpoint;

    int checkpoint_reached = (mode == OUTPUT_SCHEDULE_INTERVAL_ITERATION && iter_reached) ||
                             (mode == OUTPUT_SCHEDULE_INTERVAL_TIME && stime_reached) ||
                             (mode == OUTPUT_SCHEDULE_LIST_ITERATION && iter_reached) ||
                             (mode == OUTPUT_SCHEDULE_LIST_TIME && stime_reached);

    if (checkpoint_reached) {
        if (mode == OUTPUT_SCHEDULE_INTERVAL_ITERATION || mode == OUTPUT_SCHEDULE_INTERVAL_TIME) {
            // interval-based increments
            if (!checkpoint_reached && (ss->iter > sched->next_checkpoint)) {
                fprintf(stderr, "Iterations (%lu) exceeded log checkpoint (%lf) without noticing\n",
                        ss->iter, sched->next_checkpoint);
                clean_and_error(EXIT_FAILURE);
            }
            if (mode == OUTPUT_SCHEDULE_INTERVAL_TIME) {
                // kmc time steps can be large and skip over checkpoints, so advance in a loop until
                // we are past the current time
                while (sched->next_checkpoint <= ss->elapsed_stime) {
                    sched->next_checkpoint += sched->interval;
                }
            } else {
                sched->next_checkpoint += sched->interval;
            }
        } else if (mode == OUTPUT_SCHEDULE_LIST_ITERATION || mode == OUTPUT_SCHEDULE_LIST_TIME) {
            // list-based increments
            sched->list_idx++;
            if (sched->list_idx < sched->list_len) {
                sched->next_checkpoint = sched->list[sched->list_idx];
                if (mode == OUTPUT_SCHEDULE_INTERVAL_TIME) {
                    // kmc time steps can be large and skip over checkpoints, so advance in a loop
                    // until we are past the current time
                    while (sched->next_checkpoint <= ss->elapsed_stime) {
                        fprintf(stderr,
                                "Warning: simulation time %.4lf exceeded log checkpoint %.4lf "
                                "without logging\n",
                                ss->elapsed_stime, sched->next_checkpoint);
                        sched->list_idx++;
                        sched->next_checkpoint += sched->list[sched->list_idx];
                    }
                } else {
                    sched->next_checkpoint += sched->interval;
                }
            } else {
                // end of list reached, disable further logging
                *state = false;
            }
        }
    }
    return checkpoint_reached;
}

/**
 * @brief creates suffix for xyz file based on output schedule mode and checkpoint value (iteration
 * or time)
 *
 * @param suffix
 * @param mode
 * @param checkpoint
 */
void write_xyz_suffix(char *suffix, OutputScheduleMode mode, double checkpoint)
{
    if ((mode == OUTPUT_SCHEDULE_INTERVAL_ITERATION) || (mode == OUTPUT_SCHEDULE_LIST_ITERATION)) {
        snprintf(suffix, BUFFER_SIZE, "i%lu", (unsigned long)checkpoint);
    } else if ((mode == OUTPUT_SCHEDULE_INTERVAL_TIME) || (mode == OUTPUT_SCHEDULE_LIST_TIME)) {
        snprintf(suffix, BUFFER_SIZE, "t%.4lf", checkpoint);
    }
    return;
}

/**
 * @brief writes to log files and increments framenums
 *
 * @param output_csv whether to output to csv
 * @param output_xyz whether to output to xyz
 * @param ss
 * @param se
 * @param ls
 */
void write_logs(const StepData *step_data, struct SimulationState *ss, struct SimulationEnv *se,
                struct LoggingState *ls)
{
    // update framenums after, so initial logs (t=0) show frame 0
    int logging = 0;
    for (int i = 0; i < ls->out_formats_cnt; i++) {
        OutputFormat *format = &(ls->out_formats[i]);
        if (format->should_log_now) {
            switch (format->type) {
            case OUTPUT_FORMAT_CSV:
                logging |= format->should_log_now;
                log_state_csv(&format->csv, ls->stime_precision, ls->overpot_precision, ss);
                format->csv.frame_num++;
                break;
            case OUTPUT_FORMAT_XYZ:
                logging |= format->should_log_now;
                // suffix is expected to be updated by caller
                write_xyz_file(format->xyz.prefix, format->xyz.frame_num, format->xyz.suffix,
                               format->xyz.stripped, ss, se);
                format->xyz.frame_num++;
            case OUTPUT_FORMAT_STEPS_CSV:
                // steps csv is written to every step
                switch (se->flavor) {
                case FLAVOR_KMC:
                    log_kmc_steps(format->steps.file, step_data, ls->stime_precision);
                    break;
                case FLAVOR_MC:
                    log_mc_steps(format->steps.file, step_data);
                }
                break;
            }
        }
    }

    // if something is being logged, then output the frame to sim_log
    if (logging) {
        output_log_file(ls->sim_log, ls->framenum, ss->iter, ss->elapsed_stime, ss->temperature,
                        ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
        ls->framenum++;
    }

    return;
}

/**
 * @brief outputs to log file, csv, and xyz if checkpoints are reached; updates checkpoints for next
 * output; updates framenums for next output
 *
 * @param ss
 * @param se
 * @param ls
 */
void output_on_schedule(StepData *step_data, struct SimulationState *ss, struct SimulationEnv *se,
                        struct LoggingState *ls)
{
    for (int i = 0; i < ls->out_formats_cnt; i++) {
        OutputFormat *format = &(ls->out_formats[i]);
        int output_now = false;
        switch (format->type) {
        case OUTPUT_FORMAT_CSV:
            output_now =
                check_and_advance_checkpoint(&format->csv.schedule, &format->is_active, ss);
            break;
        case OUTPUT_FORMAT_XYZ:
            write_xyz_suffix(format->xyz.suffix, format->xyz.schedule.mode,
                             format->xyz.schedule.next_checkpoint);
            output_now =
                check_and_advance_checkpoint(&format->xyz.schedule, &format->is_active, ss);
            break;
        case OUTPUT_FORMAT_STEPS_CSV:
            // steps csv is written to every step, so no checkpoint to check against
            break;
        }
    }
    write_logs(step_data, ss, se, ls);
    return;
}
