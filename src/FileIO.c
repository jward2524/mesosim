#include "FileIO.h"
#include "Atoms.h"
#include "ErrorM.h"
#include "Input.h"
#include "Utils.h"
#include "InputXYZ.h"
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
char outFile[260] = ""; // MAX_PATH variable Windows related, default 260

static void calloc_nnE(struct SimulationEnv *se);
static int parse_datalog_params(char *params, int cursor, struct LoggingState *ls, FILE *temp_log);
static void parse_log_list(char *input_str, double *list, int *len);

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
    }

    int fres = fflush(stream);
    if (fres == EOF) {
        perror("Failed to flush log");
        clean_and_error(errno);
    }
}

static void open_log_files(struct LoggingState *ls, unsigned flavor)
{
    if (ls->output_state_csv) {
        ls->state_csv = fopen(ls->csv_filename, "w+");
        fopen_error(ls->csv_filename, ls->state_csv, "Failed to open csv file, ");

        if (flavor == FLAVOR_KMC) {
            output_kmc_csv_header(ls->state_csv);
        } else if (flavor == FLAVOR_MC) {
            output_mc_csv_header(ls->state_csv);
        }
    } else {
        ls->state_csv = NULL;
    }

    if (ls->output_steps_csv) {
        ls->iter_csv = fopen(ls->steps_filename, "w+");
        fopen_error(ls->steps_filename, ls->iter_csv, "Failed to open steps csv file, ");

        if (flavor == FLAVOR_KMC) {
            output_kmc_iter_header(ls->iter_csv);
        } else if (flavor == FLAVOR_MC) {
            output_mc_iter_header(ls->iter_csv);
        }
    } else {
        ls->iter_csv = NULL;
    }
}

bool simulation_parameters_from_file(char *filename, struct SimulationState *ss,
                                     struct SimulationEnv *se, struct LoggingState *ls)
{
    FILE *input_file = fopen(filename, "r");
    fopen_error(filename, input_file, "Failed to open input file, ");

    process_in_file(input_file, ss, se, ls);
    open_log_files(ls, se->flavor);
    fclose(input_file);
    return true;
}

void process_in_file(FILE *input_file, struct SimulationState *ss, struct SimulationEnv *se,
                     struct LoggingState *ls)
{
    ParseContext ctx = {0};
    parse_input_file(input_file, &ctx, ss, se, ls);
    finalize_config(&ctx, ss, se, ls);
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
bool process_xyz_file(FILE *input_file, struct SimulationState *ss,
                      struct SimulationEnv *se, struct LoggingState *ls)
{
    // processes file with .xyz format (number of atoms / comment / type x y z)

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
        } else if (strncmp(kv.key, "frame", 5) == 0) {
            ls->framenum = atoi(kv.value);
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

const char* get_iteration(const struct SimulationState *ss) {
    int len = snprintf(NULL, 0, "%lu", ss->iter);
    char *buf = (char*)malloc(len + 1);
    if (!buf) {
        return NULL;
    }
    snprintf(buf, len + 1, "%lu", ss->iter);
    return buf;
}

const char* get_time(const struct SimulationState *ss) {
    int len = snprintf(NULL, 0, "%le", ss->elapsed_stime);
    char *buf = (char*)malloc(len + 1);
    if (!buf) {
        return NULL;
    }
    snprintf(buf, len + 1, "%le", ss->elapsed_stime);
    return buf;
}

const char* get_energy(const struct SimulationState *ss) {
    int len = snprintf(NULL, 0, "%lf", ss->total_internal_energy);
    char *buf = (char*)malloc(len + 1);
    if (!buf) {
        return NULL;
    }
    snprintf(buf, len + 1, "%lf", ss->total_internal_energy);
    return buf;
}

const char* get_temperature(const struct SimulationState *ss) {
    int len = snprintf(NULL, 0, "%lf", ss->temperature);
    char *buf = (char*)malloc(len + 1);
    if (!buf) {
        return NULL;
    }
    snprintf(buf, len + 1, "%lf", ss->temperature);
    return buf;
}

const char* get_overpotential(const struct SimulationState *ss) {
    int len = snprintf(NULL, 0, "%lf", ss->overpotential);
    char *buf = (char*)malloc(len + 1);
    if (!buf) {
        return NULL;
    }
    snprintf(buf, len + 1, "%lf", ss->overpotential);
    return buf;
}

const char* get_total_atoms_dissolved(const struct SimulationState *ss) {
    int len = snprintf(NULL, 0, "%d", ss->total_atoms_dissolved);
    char *buf = (char*)malloc(len + 1);
    if (!buf) {
        return NULL;
    }
    snprintf(buf, len + 1, "%d", ss->total_atoms_dissolved);
    return buf;
}

const char* get_mmc_steps(const struct SimulationState *ss) {
    int len = snprintf(NULL, 0, "%lu", ss->iter);
    char *buf = (char*)malloc(len + 1);
    if (!buf) {
        return NULL;
    }
    snprintf(buf, len + 1, "%lu", ss->iter);
    return buf;
}

const CsvFieldFunc csv_field_map[] = {
    {"iter", get_iteration, FLAVOR_UNDEFINED},
    {"time", get_time, FLAVOR_UNDEFINED},
    {"energy", get_energy, FLAVOR_UNDEFINED},
    {"temperature", get_temperature, FLAVOR_UNDEFINED},
    {"overpotential", get_overpotential, FLAVOR_KMC},
    {"total_atoms_dissolved", get_total_atoms_dissolved, FLAVOR_KMC},
    {"mmc_steps", get_mmc_steps, FLAVOR_MC}
};

const size_t CSV_FIELD_FUNCS_COUNT = sizeof(csv_field_map) / sizeof(CsvFieldFunc);

// print a lot of information to the log
void input_logging(struct SimulationState *ss, struct SimulationEnv *se, struct LoggingState *ls)
{
    safe_log(ls->sim_log, "Successfully read input file and preprocessed\n");
    safe_log(ls->sim_log, "System size is %d x %d x %d\n", se->system_size_x, se->system_size_y,
             se->system_size_z);

    safe_log(ls->sim_log, "Crystal structure is ");
    switch (se->lattice_type) {
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

    safe_log(ls->sim_log, "Temperature is %lf K\n", ss->temperature);

    if (se->overpotential_ramp_rate > 0.)
        safe_log(ls->sim_log, "Potential sweep [eV/s] from %lf to %lf at %lf\n", ss->overpotential,
                 se->overpotential_ramp_rate, se->max_overpotential);
    else
        safe_log(ls->sim_log, "Potential constant [eV] at %lf\n", ss->overpotential);

    if (ls->analysis_type == REGULAR_TIME_INTERVALS)
        safe_log(ls->sim_log,
                 "Recording data at linear intervals [s] from %lf to %lf at %lf increments\n",
                 ls->next_log_checkpoint, ss->run_stime, ls->log_interval);
    else if (ls->analysis_type == LN_TIME_INTERVALS)
        safe_log(ls->sim_log,
                 "Recording data at log intervals [s] from %lf to %lf at %lf multiples\n",
                 ls->next_log_checkpoint, ss->run_stime, ls->log_interval);
    // TODO: fill out for other analysis_types

    safe_log(ls->sim_log, "Random seed is %u\n", se->rand_seed);

    switch (se->geometry) {
    case GEOMETRY_FLAT_SHEET:
        safe_log(ls->sim_log, "Initialized flat sheet with monolayer depth %d\n",
                 se->sheet_thickness);
        break;
    case GEOMETRY_CLUSTER:
        safe_log(ls->sim_log, "Initialized spherical cluster with radius %d\n", se->cluster_radius);
        break;
    case GEOMETRY_FROM_FILE:
        safe_log(ls->sim_log, "Initialized user-defined structure with filename %s\n",
                 se->atoms_filename);
        break;
    }

    safe_log(ls->sim_log, "Atoms created, %ld total\n", ss->atom_cnt);
}

/* === Data output === */

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

void output_kmc_csv_header(FILE *csv_file)
{
    safe_log(csv_file, "frame,iter,elapsed_stime,energy,temperature,overpotential,atoms\n");
}

void log_kmc_state_csv(FILE *csv_file, int frame_num, unsigned long int iter, double elapsed_stime,
                       double temperature, double overpotential, long int atom_cnt,
                       double total_internal_energy)
{
    safe_log(csv_file, "%d,", frame_num);
    safe_log(csv_file, "%lu,", iter);
    safe_log(csv_file, "%le,", elapsed_stime);
    safe_log(csv_file, "%lf,", total_internal_energy);
    safe_log(csv_file, "%lf,", temperature);
    safe_log(csv_file, "%lf,", overpotential);
    safe_log(csv_file, "%ld\n", atom_cnt);
}

void output_mc_csv_header(FILE *csv_file)
{
    safe_log(csv_file, "frame,mmc_step,iter,energy\n");
}

void log_mc_state_csv(FILE *csv_file, const int frame_num, const unsigned long int mmc_steps,
                      const unsigned long int iter, const double sys_energy)
{
    safe_log(csv_file, "%d,", frame_num);
    safe_log(csv_file, "%lu,", mmc_steps);
    safe_log(csv_file, "%lu,", iter);
    safe_log(csv_file, "%lf", sys_energy);
    safe_log(csv_file, "\n");
}

void output_kmc_iter_header(FILE *csv_file)
{
    safe_log(csv_file, "iter,sim_time,energy,u1,v1,w1,u2,v2,w2\n");
}

// output to csv:
// iteration number, simulation time, system energy (per atom?), x1, y1, z1, x2, y2, z2
// and atom ids at some point
void log_kmc_iter(FILE *csv_file, const unsigned long int iter, const double sim_time,
                  const double sys_energy, const int uvw1[3], const int uvw2[3], int is_evap)
{
    safe_log(csv_file, "%lu,", iter);
    safe_log(csv_file, "%le,", sim_time);
    safe_log(csv_file, "%lf,", sys_energy);
    safe_log(csv_file, "%d,%d,%d,", uvw1[0], uvw1[1], uvw1[2]);
    if (!is_evap) {
        safe_log(csv_file, "%d,%d,%d", uvw2[0], uvw2[1], uvw2[2]);
    } else {
        safe_log(csv_file, ",,");
    }
    safe_log(csv_file, "\n");
}

void output_mc_iter_header(FILE *csv_file)
{
    safe_log(csv_file, "iter,energy,deltaE,performed,u1,v1,w1,u2,v2,w2\n");
}

// output to csv:
// MCSS, system energy (per atom?), uvw1, uvw2
void log_mc_iter(FILE *csv_file, const unsigned long int iter, const double sys_energy,
                 const double deltaE, const int performed, const int uvw1[3], const int uvw2[3])
{
    safe_log(csv_file, "%lu,", iter);
    safe_log(csv_file, "%lf,", sys_energy);
    safe_log(csv_file, "%lf,", deltaE);
    safe_log(csv_file, "%d,", performed);
    safe_log(csv_file, "%d,%d,%d,", uvw1[0], uvw1[1], uvw1[2]);
    safe_log(csv_file, "%d,%d,%d", uvw2[0], uvw2[1], uvw2[2]);
    safe_log(csv_file, "\n");
}

bool write_xyz_file(char *xyz_filename, int frame_num, char *suffix, struct SimulationState *ss,
                    struct SimulationEnv *se)
{
    bool is_extended = 1;

    char filename_full[BUFFER_SIZE];
    sprintf(filename_full, "%s_%d_%s.xyz", xyz_filename, frame_num, suffix);
    FILE *file = fopen(filename_full, "w+");
    if (file == NULL) {
        printf("ERROR! Couldn't open output file %s\n", filename_full);
        fprintf(stderr, "Couldn't open file %s: %s\n", filename_full, strerror(errno));
        clean_and_error(errno);
    }

    /* format:
            [number of atoms]
            [comment line - exactly one line]
            [element] [x] [y] [z]
    */

    // start with number of atoms
    safe_log(file, "%ld\n", ss->atom_cnt);

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
    for (int i = 0; i < ss->atom_cnt; ++i) {
        safe_log(file, "%d %s %lf %lf %lf\n", i, se->atom_names[atoms[i]->type],
                 atoms[i]->cartesian[0], atoms[i]->cartesian[1],
                 atoms[i]->cartesian[2]); // name is now element type
    }
    // ball and stick or space filling?
    fclose(file);

#if (!defined(NDEBUG)) && defined(HAVE_FORK)
    fprintf(stderr, "Creating core dump for frame %d\n", frame_num);
    create_coredump();
#endif

    return true;
}
