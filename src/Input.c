#include "Input.h"
#include "Utils.h"
#include "ErrorM.h"
#include "FileIO.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_LINE   1024
#define MAX_TOKENS 256

// if function returns 0, function failed
// if it returns 1, handler succeeeded

/* ================= Utilities ================= */
// TODO: replace with InputXYZ.c:tokenize_line
// static int tokenize(char *line, char **argv) {
//     int argc = 0;
//     char *tok = strtok(line, " \t\r\n");
//     while (tok && argc < MAX_TOKENS) {
//         argv[argc++] = tok;
//         tok = strtok(NULL, " \t\r\n");
//     }
//     return argc;
// }
// 

// int tokenize_line(char *line, char **tokens, int maxtok);

static int parse_int(const char *s, int *out) {
    char *e;
    errno = 0;
    long v = strtol(s, &e, 10);
    if ((e == s) || (*e != '\0') || (errno != 0)) {
        return 1;
    }
    *out = (int)v;
    return 0;
}

static int parse_uint(const char *s, unsigned int *out) {
    char *e;
    errno = 0;
    unsigned long v = strtoul(s, &e, 10);
    if ((e == s) || (*e != '\0') || (errno != 0)) {
        return 1;
    }
    *out = (unsigned int)v;
    return 0;
}

static int parse_ulong(const char *s, unsigned long *out) {
    char *e;
    errno = 0;
    unsigned long v = strtoul(s, &e, 10);
    if ((e == s) || (*e != '\0') || (errno != 0)) {
        return 1;
    }
    *out = v;
    return 0;
}

static int parse_double(const char *s, double *out) {
    char *e;
    errno = 0;
    double v = strtod(s, &e);
    if ((e == s) || (*e != '\0')  || (errno != 0)) {
        return 1;
    }
    *out = v;
    return 0;
}

/* ================= Command handlers ================= */

static int cmd_atomtype(int argc, char **argv, int line, ParseContext *ctx,
                        struct SimulationState *ss, struct SimulationEnv *se,
                        struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)ls;
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
        return 1;
    }

    se->atom_names_cnt = argc - 1;
    se->atom_names = (char **)malloc((size_t)se->atom_names_cnt * sizeof(char *));
    for (int i = 0; i < se->atom_names_cnt; i++) {
        se->atom_names[i] = dup_str(argv[i+1]);
    }
    se->num_elements = argc - 1;
    se->num_bond_types = get_num_bond_types(se->num_elements);

    return 0;
}

static int cmd_composition(int argc, char **argv, int line,
                           ParseContext *ctx, struct SimulationState *ss,
                           struct SimulationEnv *se, struct LoggingState *ls) {
    (void)ss;
    (void)se;
    (void)ls;
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
        return 1;
    }

    ctx->composition_line = line;
    ctx->composition_count = argc - 1;
    ctx->composition_raw = malloc((size_t)ctx->composition_count * sizeof(double));

    for (int i = 0; i < ctx->composition_count; i++) {
        int res = parse_double(argv[i+1], &ctx->composition_raw[i]);
        if (res) {
            fprintf(stderr, "Input error - could not read composition %s\n", argv[i+1]);
            return 1;
        }
    }
    return 0;
}

static int cmd_dissolution(int argc, char **argv, int line,
                           ParseContext *ctx, struct SimulationState *ss,
                           struct SimulationEnv *se, struct LoggingState *ls) {
    (void)ss;
    (void)se;
    (void)ls;
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
        return 1;
    }

    ctx->dissolution_count = argc - 1;
    ctx->dissolution_line = line;
    ctx->dissolution_raw = malloc((size_t)ctx->dissolution_count * sizeof(int));

    for (int i = 0; i < ctx->dissolution_count; i++) {
        switch (parse_boolean(argv[i + 1])) {
        case 1:
            ctx->dissolution_raw[i] = 1;
            break;
        case 0:
            ctx->dissolution_raw[i] = 0;
            break;
        default:
            fprintf(stderr, "Input error - could not parse boolean value %s\n", argv[i + 1]);
            return 1;
            break;
        }
    }
    return 0;
}

static int cmd_nnlevels(int argc, char **argv, int line,
                        ParseContext *ctx, struct SimulationState *ss,
                        struct SimulationEnv *se, struct LoggingState *ls) {
    (void)line;
    (void)ctx;
    (void)ss;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return 1;
    }
    int res = parse_int(argv[1], &se->num_nn_levels);
    if (res) {
        fprintf(stderr, "Input error - could not parse nnlevels value %s\n", argv[1]);
        return 1;
    }
    return 0;
}

static int cmd_nne(int argc, char **argv, int line,
                   ParseContext *ctx, struct SimulationState *ss,
                   struct SimulationEnv *se, struct LoggingState *ls) {
    (void)ss;
    (void)se;
    (void)ls;
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
        return 1;
    }

    int level;
    int ret = sscanf(argv[0], "%dnne", &level);
    if ((ret != 1) || level <= 0) {
        fprintf(stderr, "Input error - expected nne command format of [level]nne, got %s\n", argv[0]);
        return 1;
    }

    int nvals = argc - 1;
    double *values = malloc((size_t)nvals * sizeof(double));
    for (int i = 0; i < nvals; i++) {
        int pres = parse_double(argv[i + 1], &values[i]);
        if (pres) {
            fprintf(stderr, "Input error - could not parse nne value %s\n", argv[i + 1]);
            free(values);
            return 1;
        }
    }

    ctx->nne_cmds = realloc(ctx->nne_cmds, (size_t)(ctx->nne_cmd_count + 1) * sizeof(*ctx->nne_cmds));

    ctx->nne_cmds[ctx->nne_cmd_count] = (struct NNECmd) {
        .level = level,
        .line = line,
        .values = values,
        .nvalues = nvals
    };
    ctx->nne_cmd_count++;

    return 0;
}

static int cmd_systemsize(int argc, char **argv, int line,
                          ParseContext *ctx, struct SimulationState *ss,
                          struct SimulationEnv *se, struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)ls;
    if (argc != 4) {
        fprintf(stderr, "Input error - expected 3 arguments for command %s\n", argv[0]);
        return 1;
    }
    if (parse_int(argv[1], &se->system_size_x) || parse_int(argv[2], &se->system_size_y) ||
        parse_int(argv[3], &se->system_size_z)) {
        fprintf(stderr, "Input error - could not parse systemsize arguments\n");
        return 1;
    }
    se->max_atoms = (long)se->system_size_x * (long)se->system_size_y * (long)se->system_size_z;
    return 0;
}

static int cmd_temp(int argc, char **argv, int line, ParseContext *ctx,
                    struct SimulationState *ss, struct SimulationEnv *se,
                    struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)se;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return 1;
    }
    if (parse_double(argv[1], &ss->temperature)) {
        fprintf(stderr, "Input error - could not parse temp value %s\n", argv[1]);
        return 1;
    }
    return 0;
}

static int cmd_seed(int argc, char **argv, int line, ParseContext *ctx,
                    struct SimulationState *ss, struct SimulationEnv *se,
                    struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "random") == 0) {
        time_t seedtime;
        time(&seedtime);
        se->rand_seed = (unsigned int)seedtime;
        return 0;
    }
    if (strcmp(argv[1], "default") == 0) {
        se->rand_seed = DEFAULT_SEED;
        return 0;
    }
    if (parse_uint(argv[1], &se->rand_seed)) {
        fprintf(stderr, "Input error - could not parse seed value %s\n", argv[1]);
        return 1;
    }
    return 0;
}

static int cmd_potential(int argc, char **argv, int line, ParseContext *ctx,
                         struct SimulationState *ss, struct SimulationEnv *se,
                         struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)ls;
    if (argc != 2 && argc != 4) {
        fprintf(stderr, "Input error - expected 1 or 3 arguments for command %s\n", argv[0]);
        return 1;
    }
    if (parse_double(argv[1], &se->initial_overpotential)) {
        fprintf(stderr, "Input error - could not parse potential value %s\n", argv[1]);
        return 1;
    }
    if (argc == 2) {
        se->overpotential_ramp_rate = 0.0;
        se->max_overpotential = se->initial_overpotential;
        return 0;
    }
    if (parse_double(argv[2], &se->overpotential_ramp_rate) ||
        parse_double(argv[3], &se->max_overpotential)) {
        fprintf(stderr, "Input error - could not parse potential sweep arguments\n");
        return 1;
    }
    return 0;
}

static int cmd_datalog(int argc, char **argv, int line, ParseContext *ctx,
                       struct SimulationState *ss, struct SimulationEnv *se,
                       struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)se;
    if (argc < 4) {
        fprintf(stderr, "Input error - not enough arguments for command %s\n", argv[0]);
        return 1;
    }

    int interval_type = 0;
    if (strcmp(argv[1], "lineart") == 0) {
        interval_type = REGULAR_TIME_INTERVALS;
    } else if (strcmp(argv[1], "lnt") == 0) {
        interval_type = LN_TIME_INTERVALS;
    } else if (strcmp(argv[1], "iteration") == 0) {
        interval_type = ITERATION_INTERVALS;
    } else {
        fprintf(stderr, "Input error - unknown datalog type %s\n", argv[1]);
        return 1;
    }

    if (strcmp(argv[2], "interval") == 0) {
        if (argc != 4 && argc != 5) {
            fprintf(stderr, "Input error - datalog interval expects 1 or 2 numeric values\n");
            return 1;
        }
        if (parse_double(argv[3], &ls->next_log_checkpoint)) {
            fprintf(stderr, "Input error - invalid datalog checkpoint %s\n", argv[3]);
            return 1;
        }
        if (argc == 5) {
            if (parse_double(argv[4], &ls->log_interval)) {
                fprintf(stderr, "Input error - invalid datalog interval %s\n", argv[4]);
                return 1;
            }
        } else {
            ls->log_interval = ls->next_log_checkpoint;
        }
        ls->analysis_type = interval_type;
    } else if (strcmp(argv[2], "list") == 0) {
        if (argc < 4) {
            fprintf(stderr, "Input error - datalog list requires at least one value\n");
            return 1;
        }
        int nvals = argc - 3;
        ls->log_list = (double *)malloc((size_t)nvals * sizeof(double));
        if (!ls->log_list) {
            fprintf(stderr, "Input error - failed to allocate datalog list\n");
            return 1;
        }
        for (int i = 0; i < nvals; i++) {
            if (parse_double(argv[i + 3], &ls->log_list[i])) {
                fprintf(stderr, "Input error - invalid datalog list value %s\n", argv[i + 3]);
                free(ls->log_list);
                ls->log_list = NULL;
                ls->log_list_len = 0;
                return 1;
            }
        }
        ls->log_list_len = nvals;
        ls->analysis_type = (interval_type == ITERATION_INTERVALS) ? ITERATION_LIST : TIME_LIST;
    } else {
        fprintf(stderr, "Input error - expected datalog mode interval|list, got %s\n", argv[2]);
        return 1;
    }

    ls->framenum = 0;
    return 0;
}

static int cmd_struct(int argc, char **argv, int line, ParseContext *ctx,
                      struct SimulationState *ss, struct SimulationEnv *se,
                      struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "FCC") == 0 || strcmp(argv[1], "fcc") == 0) {
        se->lattice_type = FCC;
    } else if (strcmp(argv[1], "BCC") == 0 || strcmp(argv[1], "bcc") == 0) {
        se->lattice_type = BCC;
    } else if (strcmp(argv[1], "SC") == 0 || strcmp(argv[1], "sc") == 0) {
        se->lattice_type = SC;
    } else {
        fprintf(stderr, "Input error - structure type %s not valid\n", argv[1]);
        return 1;
    }
    se->num_transition_vectors = MAXIMUM_NUMBER_OF_NEIGHBORS;
    return 0;
}

static int cmd_output(int argc, char **argv, int line, ParseContext *ctx,
                      struct SimulationState *ss, struct SimulationEnv *se,
                      struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)se;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return 1;
    }
    snprintf(outFile, 260, "%s", argv[1]);
    return 0;
}

static int cmd_geometry(int argc, char **argv, int line, ParseContext *ctx,
                        struct SimulationState *ss, struct SimulationEnv *se,
                        struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)ls;
    if (argc < 3) {
        fprintf(stderr, "Input error - not enough arguments for command %s\n", argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "sheet") == 0) {
        if (argc != 3 || parse_int(argv[2], &se->sheet_thickness)) {
            fprintf(stderr, "Input error - invalid sheet thickness %s\n", argc > 2 ? argv[2] : "");
            return 1;
        }
        se->geometry = GEOMETRY_FLAT_SHEET;
    } else if (strcmp(argv[1], "cluster") == 0) {
        if (argc != 3 || parse_int(argv[2], &se->cluster_radius)) {
            fprintf(stderr, "Input error - invalid cluster radius %s\n", argc > 2 ? argv[2] : "");
            return 1;
        }
        se->geometry = GEOMETRY_CLUSTER;
    } else if (strcmp(argv[1], "file") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Input error - geometry file expects a filename\n");
            return 1;
        }
        se->geometry = GEOMETRY_FROM_FILE;
        snprintf(se->atoms_filename, sizeof(se->atoms_filename), "%s", argv[2]);
    } else {
        fprintf(stderr, "Input error - unrecognized geometry type %s\n", argv[1]);
        return 1;
    }
    return 0;
}

static int cmd_run(int argc, char **argv, int line, ParseContext *ctx,
                   struct SimulationState *ss, struct SimulationEnv *se,
                   struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)se;
    (void)ls;
    if (argc != 3) {
        fprintf(stderr, "Input error - expected 2 arguments for command %s\n", argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "time") == 0) {
        ss->sim_end_type = SIM_END_BY_STIME;
        if (parse_double(argv[2], &ss->run_stime)) {
            fprintf(stderr, "Input error - invalid run time %s\n", argv[2]);
            return 1;
        }
    } else if (strcmp(argv[1], "iteration") == 0) {
        ss->sim_end_type = SIM_END_BY_ITERATIONS;
        if (parse_ulong(argv[2], &ss->final_iteration)) {
            fprintf(stderr, "Input error - invalid run iteration %s\n", argv[2]);
            return 1;
        }
    } else {
        fprintf(stderr, "Input error - unknown run mode %s\n", argv[1]);
        return 1;
    }
    return 0;
}

static int cmd_flavor(int argc, char **argv, int line, ParseContext *ctx,
                      struct SimulationState *ss, struct SimulationEnv *se,
                      struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "KMC") == 0) {
        se->flavor = FLAVOR_KMC;
    } else if (strcmp(argv[1], "MC") == 0) {
        se->flavor = FLAVOR_MC;
    } else {
        fprintf(stderr, "Input error - unknown flavor %s\n", argv[1]);
        return 1;
    }
    return 0;
}

static int cmd_logtype(int argc, char **argv, int line, ParseContext *ctx,
                       struct SimulationState *ss, struct SimulationEnv *se,
                       struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)se;
    if (argc < 2) {
        fprintf(stderr, "Input error - expected at least 1 argument for command %s\n", argv[0]);
        return 1;
    }
    for (int i = 1; i < argc; i++) {
        if (strstr(argv[i], "iter") != NULL) {
            ls->output_iter_csv = 1;
        }
        if (strstr(argv[i], "csv") != NULL) {
            ls->output_state_csv = 1;
        }
        if (strstr(argv[i], "xyz") != NULL) {
            ls->output_xyz = 1;
        }
    }
    return 0;
}

/* ================= Command table ================= */
static int cmd_help(int argc, char **argv, int line, ParseContext *ctx,
                    struct SimulationState *ss, struct SimulationEnv *se,
                    struct LoggingState *ls);

Command commands[] = {
    { "systemsize",  cmd_systemsize,  "systemsize NX NY NZ", "Simulation box size in lattice units.",  CMDCAT_GEOMETRY, 1},
    { "temp",        cmd_temp,        "temp T", "Simulation temperature in K.", CMDCAT_THERMODYNAMICS, 1},
    { "seed",        cmd_seed,        "seed random|default|N", "Selection of the random seed value. Default: default", CMDCAT_RUN, 0},
    { "potential",   cmd_potential,   "potential U0 [dUdt Umax]", "Constant or swept electric potential. Default: 0", CMDCAT_THERMODYNAMICS, 0},
    { "datalog",     cmd_datalog,     "datalog (linear|ln|iteration) (interval a [b]|list ...)", "Logging schedule configuration.", CMDCAT_OUTPUT, 0},
    { "struct",      cmd_struct,      "struct FCC|BCC|SC", "Crystal structure type.", CMDCAT_GEOMETRY, 1},
    { "output",      cmd_output,      "output path/outfile.out", "Output log filename.", CMDCAT_OUTPUT, 0},
    { "geometry",    cmd_geometry,    "geometry (sheet N|cluster R|file path)", "Initial geometry configuration.", CMDCAT_GEOMETRY, 1},
    { "atomtype",    cmd_atomtype,    "atomtype A B [C ...]",  "Define atom types and their order.", CMDCAT_GEOMETRY, 1},
    { "composition", cmd_composition, "composition xA xB [xC ...]", "Atomic composition fractions; order follows atomtype.", CMDCAT_GEOMETRY, 1},
    { "dissolution", cmd_dissolution, "dissolution true|false ...", "Dissolution flags per atom type; order follows atomtype.", CMDCAT_GEOMETRY, 1},
    { "nnlevels",    cmd_nnlevels,    "nnlevels N", "Number of nearest-neighbor shells.", CMDCAT_THERMODYNAMICS, 1},
    { "nne",         cmd_nne,         "1nne eAA eAB ...", "Nearest-neighbor energies for shell n (flattened upper-triangle). For a three-component system: AA AB AC BB BC CC", CMDCAT_THERMODYNAMICS, 1},
    { "run",         cmd_run,         "run time|iteration value", "Simulation end condition.", CMDCAT_RUN, 1},
    { "flavor",      cmd_flavor,      "flavor KMC|MC", "Simulation algorithm flavor.", CMDCAT_RUN, 1},
    { "logtype",     cmd_logtype,     "logtype [iter] [csv] [xyz]", "Enable output formats.", CMDCAT_OUTPUT, 0},
    { "help",        cmd_help,        "help [command]", "Show documentation for commands.", CMDCAT_OUTPUT, 0},
    { NULL, NULL, NULL, NULL }
};

/* ================= Parser ================= */

void parse_input_file(FILE *fp, ParseContext *ctx, struct SimulationState *ss,
                      struct SimulationEnv *se, struct LoggingState *ls) {
    char line[MAX_LINE];
    char *argv[MAX_TOKENS];
    int lineno = 0;

    while (fgets(line, sizeof(line), fp)) {
        lineno++;
        if (line[0] == '#')
            continue;

        int argc = tokenize_line(line, argv, MAX_TOKENS);
        if (argc == 0) continue;

        int matched = 0;
        for (Command *c = commands; c->name != NULL; c++) {
            if (strcmp(argv[0], c->name) == 0 ||
                (strcmp(c->name, "nne") == 0 && isdigit(argv[0][0]))) {
                int ret = c->func(argc, argv, lineno, ctx, ss, se, ls);
                if (ret) {
                    // if handler fails
                    fprintf(stderr,
                        "Error in command '%s' at line %d\nUsage: %s\n",
                        argv[0], lineno, c->usage);
                    call_exit(EXIT_FAILURE);
                }
                matched = 1;
                break;
            }
        }

        if (!matched) {
            fprintf(stderr, "Unknown command '%s' at line %d\n",
                    argv[0], lineno);
            call_exit(EXIT_FAILURE);
        }
    }
}

/* ================= Help ================= */

void print_help(const char *cmd) {
    for (Command *c = commands; c->name != NULL; c++) {
        if (!cmd || strcmp(cmd, c->name) == 0) {
            printf("%s\n  %s\n\n", c->usage, c->description);
        }
    }
}

static int cmd_help(int argc, char **argv, int line, ParseContext *ctx,
                    struct SimulationState *ss, struct SimulationEnv *se,
                    struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ss;
    (void)se;
    (void)ls;
    if (argc == 1)
        print_help(NULL);
    else
        print_help(argv[1]);
    return 1;
}


/* ================= Finalization ================= */

void finalize_atom_dependent(const ParseContext *ctx, struct SimulationEnv *se)
{
    if (se->atom_names_cnt <= 0) {
        fprintf(stderr, "atomtype must be specified\n");
        clean_and_error(EXIT_FAILURE);
    }

    if (ctx->composition_raw) {
        if (ctx->composition_count != se->atom_names_cnt) {
            fprintf(stderr, "composition count mismatch at line %d\n", ctx->composition_line);
            clean_and_error(EXIT_FAILURE);
        }
        double total = 0.0;
        for (int i = 0; i < ctx->composition_count; i++) {
            total += ctx->composition_raw[i];
        }
        if (fabs(total - 1.0) > 1e-10) {
            fprintf(stderr, "composition values must sum to 1.0, got %.17g\n", total);
            clean_and_error(EXIT_FAILURE);
        }
        se->substrate_composition = (double *)malloc((size_t)se->atom_names_cnt * sizeof(double));
        memcpy(se->substrate_composition, ctx->composition_raw,
               (size_t)se->atom_names_cnt * sizeof(double));
    }

    if (ctx->dissolution_raw) {
        if (ctx->dissolution_count != se->atom_names_cnt) {
            fprintf(stderr, "dissolution count mismatch at line %d\n", ctx->dissolution_line);
            clean_and_error(EXIT_FAILURE);
        }
        se->is_soluble = (bool *)malloc((size_t)se->atom_names_cnt * sizeof(bool));
        memcpy(se->is_soluble, ctx->dissolution_raw, (size_t)se->atom_names_cnt * sizeof(bool));
        se->dissolution = 0;
        for (int i = 0; i < ctx->dissolution_count; i++) {
            if (ctx->dissolution_raw[i]) {
                se->dissolution = 1;
                break;
            }
        }
    }
}

void finalize_nne(const ParseContext *ctx, struct SimulationEnv *se)
{
    int ntypes = se->atom_names_cnt;
    se->num_bond_types = get_num_bond_types(se->num_elements);

    if (se->num_nn_levels <= 0) {
        fprintf(stderr, "nnlevels must be specified\n");
        clean_and_error(EXIT_FAILURE);
    }

    // Allocate 3D array nne[shell][i][j]
    // se->nn_energy = malloc(se->num_nn_levels * sizeof(double**));
    // for (int s = 0; s < se->num_nn_levels; s++) {
    //     se->nn_energy[s] = malloc(ntypes * sizeof(double*));
    //     for (int i = 0; i < ntypes; i++)
    //         se->nn_energy[s][i] = calloc(ntypes, sizeof(double));
    // }

    se->num_nn_types = se->num_nn_levels * se->num_bond_types;
    se->nn_energy = (double *)calloc((size_t)se->num_nn_types, sizeof(double));
    
    int *defined = calloc((size_t)se->num_nn_levels, sizeof(int));
    int expected = se->num_bond_types;
    
    for (int c = 0; c < ctx->nne_cmd_count; c++) {
        int lvl = ctx->nne_cmds[c].level;
        if (ctx->nne_cmds[c].nvalues != expected) {
            fprintf(stderr, "nne at line %d expects %d values, got %d\n", ctx->nne_cmds[c].line,
                    expected, ctx->nne_cmds[c].nvalues);
            clean_and_error(EXIT_FAILURE);
        }
        if (lvl > se->num_nn_levels) {
            fprintf(stderr, "nne level %d exceeds nnlevels (line %d)\n", lvl,
                    ctx->nne_cmds[c].line);
            clean_and_error(EXIT_FAILURE);
        }
        if (defined[lvl-1]) {
            fprintf(stderr, "duplicate nne for level %d (line %d)\n", lvl, ctx->nne_cmds[c].line);
            clean_and_error(EXIT_FAILURE);
        }

        if (ctx->nne_cmds[c].values == NULL) {
            fprintf(stderr, "nne command at line %d has no values\n", ctx->nne_cmds[c].line);
            clean_and_error(EXIT_FAILURE);
        }

        // Fill upper triangle from flattened input
        int bond_index = nn_bondidx_2_envidx(lvl - 1, 0, se->num_bond_types);
        double *vals = ctx->nne_cmds[c].values;
        int idx = 0;
        for (int i = 0; i < ntypes; i++) {
            for (int j = i; j < ntypes; j++) {
                if (idx >= expected) break;
                se->nn_energy[bond_index] = vals[idx];
                idx++;
                bond_index++;
            }
        }
        defined[lvl-1] = 1;
    }

    for (int i = 0; i < se->num_nn_levels; i++) {
        if (!defined[i]) {
            fprintf(stderr,
                "missing nne definition for level %d\n", i+1);
            clean_and_error(EXIT_FAILURE);
        }
    }

    free(defined);
}

// static void finalize_nne(const ParseContext *ctx, struct SimulationState *ss,
//                          struct SimulationEnv *se, struct LoggingState *ls) {
//     int ntypes = se->atom_names_cnt;

//     if (se->num_nn_levels <= 0) {
//         fprintf(stderr, "nnlevels must be specified\n");
//         clean_and_error(EXIT_FAILURE);
//     }

//     // Allocate 3D array nne[shell][i][j]
//     se->nn_energy = malloc(se->num_nn_levels * sizeof(double**));
//     for (int s = 0; s < se->num_nn_levels; s++) {
//         se->nn_energy[s] = malloc(ntypes * sizeof(double*));
//         for (int i = 0; i < ntypes; i++)
//             se->nn_energy[s][i] = calloc(ntypes, sizeof(double));
//     }

//     int *defined = calloc(se->num_nn_levels, sizeof(int));

//     int expected = ntypes * (ntypes + 1) / 2;
    
//     for (int c = 0; c < ctx->nne_cmd_count; c++) {
//         int lvl = ctx->nne_cmds[c].level;
//         if (ctx->nne_cmds[c].values_count != expected) {
//             fprintf(stderr,
//                 "nne at line %d expects %d values, got %d\n",
//                 ctx->nne_cmds[c].line, expected,
//                 ctx->nne_cmds[c].values_count);
//             exit(EXIT_FAILURE);
//         }
//         if (lvl > sim->nnlevels) {
//             fprintf(stderr,
//                 "nne level %d exceeds nnlevels (line %d)\n",
//                 lvl, ctx->nne_cmds[c].line);
//             exit(EXIT_FAILURE);
//         }
//         if (defined[lvl-1]) {
//             fprintf(stderr,
//                 "duplicate nne for level %d (line %d)\n",
//                 lvl, ctx->nne_cmds[c].line);
//             exit(EXIT_FAILURE);
//         }

//         // Compute expected number of values
//         int expected = ntypes*(ntypes+1)/2;
//         if (ctx->nne_cmds[c].values == NULL) {
//             fprintf(stderr,
//                 "nne command at line %d has no values\n",
//                 ctx->nne_cmds[c].line);
//             exit(EXIT_FAILURE);
//         }

//         // Fill upper triangle from flattened input
//         double *vals = ctx->nne_cmds[c].values;
//         int idx = 0;
//         for (int i = 0; i < ntypes; i++) {
//             for (int j = i; j < ntypes; j++) {
//                 if (idx >= expected) break;
//                 sim->nne[lvl-1][i][j] = vals[idx];
//                 sim->nne[lvl-1][j][i] = vals[idx]; // symmetry
//                 idx++;
//             }
//         }
//         defined[lvl-1] = 1;
//     }

//     for (int i = 0; i < sim->nnlevels; i++) {
//         if (!defined[i]) {
//             fprintf(stderr,
//                 "missing nne definition for level %d\n", i+1);
//             exit(EXIT_FAILURE);
//         }
//     }

//     free(defined);
// }

void finalize_config(const ParseContext *ctx, struct SimulationState *ss,
                     struct SimulationEnv *se, struct LoggingState *ls) {
    (void)ss;
    (void)ls;
    finalize_atom_dependent(ctx, se);
    finalize_nne(ctx, se);
}
