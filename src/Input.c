#include "Input.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "Utils.h"
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_LINE 1024
#define MAX_TOKENS 256

// if function returns 0, function failed
// if it returns 1, handler succeeeded

/* ================= Utilities ================= */

static int parse_int(const char *s, int *out)
{
    char *e;
    errno = 0;
    long v = strtol(s, &e, 10);
    if ((e == s) || (*e != '\0') || (errno != 0)) {
        return 1;
    }
    *out = (int)v;
    return 0;
}

static int parse_uint(const char *s, unsigned int *out)
{
    char *e;
    errno = 0;
    unsigned long v = strtoul(s, &e, 10);
    if ((e == s) || (*e != '\0') || (errno != 0)) {
        return 1;
    }
    *out = (unsigned int)v;
    return 0;
}

static int parse_ulong(const char *s, unsigned long *out)
{
    char *e;
    errno = 0;
    unsigned long v = strtoul(s, &e, 10);
    if ((e == s) || (*e != '\0') || (errno != 0)) {
        return 1;
    }
    *out = v;
    return 0;
}

static int parse_double(const char *s, double *out)
{
    char *e;
    errno = 0;
    double v = strtod(s, &e);
    if ((e == s) || (*e != '\0') || (errno != 0)) {
        return 1;
    }
    *out = v;
    return 0;
}

/* ================= Command handlers ================= */

static InputErrorFlag cmd_atomtype(int argc, char **argv, int line, ParseContext *p_ctx,
                                   struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }

    inputs->atom_names_cnt = argc - 1;
    inputs->atom_names = (char **)malloc((size_t)inputs->atom_names_cnt * sizeof(char *));
    if (inputs->atom_names == NULL) {
        fprintf(stderr, "Couldn't allocate memory for atom names: %s", strerror(errno));
        return INPUT_ERR_ALLOC;
    }
    for (int i = 0; i < inputs->atom_names_cnt; i++) {
        inputs->atom_names[i] = dup_str(argv[i + 1]);
    }
    inputs->num_elements = argc - 1;

    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_composition(int argc, char **argv, int line, ParseContext *p_ctx,
                                      struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)inputs;
    (void)ls;
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }

    p_ctx->composition_lineno = line;
    p_ctx->composition_count = argc - 1;
    p_ctx->composition_raw = malloc((size_t)p_ctx->composition_count * sizeof(double));
    if (p_ctx->composition_raw == NULL) {
        fprintf(stderr, "Couldn't allocate memory for composition array in ParseContext: %s",
                strerror(errno));
        return INPUT_ERR_ALLOC;
    }

    for (int i = 0; i < p_ctx->composition_count; i++) {
        int res = parse_double(argv[i + 1], &p_ctx->composition_raw[i]);
        if (res) {
            fprintf(stderr, "Input error - could not read composition %s\n", argv[i + 1]);
            return INPUT_ERR_INVALID_ARG;
        }
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_dissolution(int argc, char **argv, int line, ParseContext *p_ctx,
                                      struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)inputs;
    (void)ls;
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }

    p_ctx->dissolution_count = argc - 1;
    p_ctx->dissolution_lineno = line;
    p_ctx->dissolution_raw = malloc((size_t)p_ctx->dissolution_count * sizeof(int));
    if (p_ctx->dissolution_raw == NULL) {
        fprintf(stderr, "Couldn't allocate memory for dissolution array in ParseContext: %s",
                strerror(errno));
        return INPUT_ERR_ALLOC;
    }

    for (int i = 0; i < p_ctx->dissolution_count; i++) {
        switch (parse_boolean(argv[i + 1])) {
        case 1:
            p_ctx->dissolution_raw[i] = 1;
            break;
        case 0:
            p_ctx->dissolution_raw[i] = 0;
            break;
        default:
            fprintf(stderr, "Input error - could not parse boolean value %s\n", argv[i + 1]);
            return INPUT_ERR_INVALID_ARG;
        }
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_nnlevels(int argc, char **argv, int line, ParseContext *p_ctx,
                                   struct SimulationConfig *inputs, struct LoggingState *ls)
{
    // TODO: this is redundant since nne command can track levels
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }
    int res = parse_int(argv[1], &inputs->num_nn_levels);
    if (res) {
        fprintf(stderr, "Input error - could not parse nnlevels value %s\n", argv[1]);
        return INPUT_ERR_INVALID_ARG;
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_nne(int argc, char **argv, int line, ParseContext *p_ctx,
                              struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)inputs;
    (void)ls;
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }

    int level;
    int ret = sscanf(argv[0], "%dnne", &level);
    if ((ret != 1) || level <= 0) {
        fprintf(stderr, "Input error - expected nne command format of [level]nne, got %s\n",
                argv[0]);
        return INPUT_ERR_INVALID_ARG;
    }

    int nvals = argc - 1;
    double *values = malloc((size_t)nvals * sizeof(double));
    if (values == NULL) {
        fprintf(stderr, "Couldn't allocate memory for nne values: %s", strerror(errno));
        return INPUT_ERR_ALLOC;
    }
    for (int i = 0; i < nvals; i++) {
        int pres = parse_double(argv[i + 1], &values[i]);
        if (pres) {
            fprintf(stderr, "Input error - could not parse nne value %s\n", argv[i + 1]);
            free(values);
            return INPUT_ERR_INVALID_ARG;
        }
    }

    p_ctx->nne_cmds =
        realloc(p_ctx->nne_cmds, (size_t)(p_ctx->nne_cmd_count + 1) * sizeof(*p_ctx->nne_cmds));
    if (p_ctx->nne_cmds == NULL) {
        fprintf(stderr, "Couldn't allocate memory for nne commands in ParseContext: %s",
                strerror(errno));
        free(values);
        return INPUT_ERR_ALLOC;
    }

    p_ctx->nne_cmds[p_ctx->nne_cmd_count] =
        (struct NNECmd){.level = level, .line = line, .values = values, .nvalues = nvals};
    p_ctx->nne_cmd_count++;

    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_systemsize(int argc, char **argv, int line, ParseContext *p_ctx,
                                     struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc != 4) {
        fprintf(stderr, "Input error - expected 3 arguments for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }
    if (parse_int(argv[1], &inputs->system_size_x) || parse_int(argv[2], &inputs->system_size_y) ||
        parse_int(argv[3], &inputs->system_size_z)) {
        fprintf(stderr, "Input error - could not parse systemsize arguments\n");
        return INPUT_ERR_INVALID_ARG;
    }
    // TODO: move to finalization
    // se->max_atoms = (long)se->system_size_x * (long)se->system_size_y * (long)se->system_size_z;
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_temp(int argc, char **argv, int line, ParseContext *ctx,
                               struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)ctx;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }
    if (parse_double(argv[1], &inputs->temperature)) {
        fprintf(stderr, "Input error - could not parse temp value %s\n", argv[1]);
        return INPUT_ERR_INVALID_ARG;
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_seed(int argc, char **argv, int line, ParseContext *p_ctx,
                               struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }

    if (strcmp(argv[1], "random") == 0) {
        time_t seedtime;
        time(&seedtime);
        inputs->rand_seed = (unsigned int)seedtime;
        return INPUT_ERR_NONE;
    }
    if (strcmp(argv[1], "default") == 0) {
        inputs->rand_seed = DEFAULT_SEED;
        return INPUT_ERR_NONE;
    }
    if (parse_uint(argv[1], &inputs->rand_seed)) {
        fprintf(stderr, "Input error - could not parse seed value %s\n", argv[1]);
        return INPUT_ERR_INVALID_ARG;
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_potential(int argc, char **argv, int line, ParseContext *p_ctx,
                                    struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc != 2 && argc != 4) {
        fprintf(stderr, "Input error - expected 1 or 3 arguments for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }
    if (parse_double(argv[1], &inputs->initial_overpotential)) {
        fprintf(stderr, "Input error - could not parse potential value %s\n", argv[1]);
        return INPUT_ERR_INVALID_ARG;
    }
    if (argc == 2) {
        inputs->overpotential_ramp_rate = 0.0;
        inputs->max_overpotential = inputs->initial_overpotential;
        return INPUT_ERR_NONE;
    }
    if (parse_double(argv[2], &inputs->overpotential_ramp_rate) ||
        parse_double(argv[3], &inputs->max_overpotential)) {
        fprintf(stderr, "Input error - could not parse potential sweep arguments\n");
        return INPUT_ERR_INVALID_ARG;
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_struct(int argc, char **argv, int line, ParseContext *p_ctx,
                                 struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }
    if (strcmp(argv[1], "FCC") == 0 || strcmp(argv[1], "fcc") == 0) {
        inputs->lattice_type = FCC;
    } else if (strcmp(argv[1], "BCC") == 0 || strcmp(argv[1], "bcc") == 0) {
        inputs->lattice_type = BCC;
    } else if (strcmp(argv[1], "SC") == 0 || strcmp(argv[1], "sc") == 0) {
        inputs->lattice_type = SC;
    } else {
        fprintf(stderr, "Input error - structure type %s not valid\n", argv[1]);
        return INPUT_ERR_INVALID_ARG;
    }
    // TODO: move to finalization
    // inputs->num_transition_vectors = MAXIMUM_NUMBER_OF_NEIGHBORS;
    return INPUT_ERR_NONE;
}

/* === Helpers for cmd_output === */

static InputErrorFlag parse_output_schedule(int argc, char **argv, int *pidx, OutputSchedule *sched,
                                            int expect_fields)
{
    // TODO: forbid time-based schedules for MC simulations
    int idx = *pidx;
    if (idx >= argc) {
        fprintf(stderr, "Input error - output expects schedule mode\n");
        return INPUT_ERR_COUNT_MISMATCH;
    }
    if (strcmp(argv[idx], "interval") == 0) {
        idx++;
        if (idx + 1 >= argc) {
            fprintf(stderr, "Input error - output interval expects type and value\n");
            return INPUT_ERR_COUNT_MISMATCH;
        }
        if (strcmp(argv[idx], "iteration") == 0) {
            sched->mode = OUTPUT_SCHEDULE_INTERVAL_ITERATION;
        } else if (strcmp(argv[idx], "time") == 0) {
            sched->mode = OUTPUT_SCHEDULE_INTERVAL_TIME;
        } else {
            fprintf(stderr, "Input error - output interval expects 'iteration' or 'time'\n");
            return INPUT_ERR_INVALID_ARG;
        }
        idx++;
        char *endptr = NULL;
        double val = strtod(argv[idx], &endptr);
        if (!endptr || *endptr != '\0') {
            fprintf(stderr, "Input error - output interval value not numeric\n");
            return INPUT_ERR_INVALID_ARG;
        }
        sched->interval = val;
        sched->list = NULL;
        sched->list_len = 0;
        sched->next_checkpoint = val;
        idx++;
    } else if (strcmp(argv[idx], "list") == 0) {
        idx++;
        if (idx >= argc) {
            fprintf(stderr, "Input error - output list expects type\n");
            return INPUT_ERR_COUNT_MISMATCH;
        }
        if (strcmp(argv[idx], "iteration") == 0) {
            sched->mode = OUTPUT_SCHEDULE_LIST_ITERATION;
        } else if (strcmp(argv[idx], "time") == 0) {
            sched->mode = OUTPUT_SCHEDULE_LIST_TIME;
        } else {
            fprintf(stderr, "Input error - output list expects 'iteration' or 'time'\n");
            return INPUT_ERR_INVALID_ARG;
        }
        idx++;

        // Parse list of values until 'fields' keyword (if expect_fields, aka if csv) or end (if
        // not, aka xyz)
        int list_start = idx;
        int list_count = 0;
        while (idx < argc && (!expect_fields || strcmp(argv[idx], "fields") != 0)) {
            list_count++;
            idx++;
        }
        if (list_count == 0) {
            fprintf(stderr, "Input error - output list expects at least one value\n");
            return INPUT_ERR_COUNT_MISMATCH;
        }

        sched->list = (double *)malloc(sizeof(double) * (size_t)list_count);
        if (!sched->list) {
            fprintf(stderr, "Couldn't allocate memory for schedule list: %s", strerror(errno));
            return INPUT_ERR_ALLOC;
        }

        for (int i = 0; i < list_count; ++i) {
            char *endptr = NULL;
            double val = strtod(argv[list_start + i], &endptr);
            if (!endptr || *endptr != '\0') {
                free(sched->list);
                sched->list = NULL;
                fprintf(stderr, "Input error - output list value not numeric\n");
                return INPUT_ERR_INVALID_ARG;
            }
            sched->list[i] = val;
        }
        sched->list_len = list_count;
        sched->list_idx = 0;
        sched->interval = 0;
        sched->next_checkpoint = sched->list[0];
        idx = list_start + list_count;
    } else {
        fprintf(stderr, "Input error - output expects 'interval' or 'list'\n");
        return INPUT_ERR_INVALID_ARG;
    }
    *pidx = idx;
    return INPUT_ERR_NONE;
}

static InputErrorFlag parse_output_csv(int argc, char **argv, struct LoggingState *ls)
{
    int idx = 2;

    OutputFormat *format = &(ls->out_formats[ls->out_formats_cnt - 1]);
    format->type = OUTPUT_FORMAT_CSV;

    // if filename is present and not a schedule keyword, use it; else generate default
    if (idx < argc && strcmp(argv[idx], "interval") != 0 && strcmp(argv[idx], "list") != 0) {
        snprintf(format->csv.filename, sizeof(format->csv.filename), "%s", argv[idx]);
        idx++;
    } else {
        // generate default filename: [time in seconds].csv
        time_t now = time(NULL);
        snprintf(format->csv.filename, sizeof(format->csv.filename), "%ld.csv", (long)now);
    }
    InputErrorFlag err = parse_output_schedule(argc, argv, &idx, &format->csv.schedule, 1);
    if (err != INPUT_ERR_NONE) {
        return err;
    }

    // parse 'fields' keyword
    // could store in ctx before finalizing, but not really necessary since fails fast
    if (idx >= argc || strcmp(argv[idx], "fields") != 0) {
        fprintf(stderr, "Input error - output csv expects 'fields' keyword\n");
        return INPUT_ERR_MISSING_CMD;
    }
    idx++;
    int field_count = argc - idx;
    if (field_count <= 0) {
        fprintf(stderr, "Input error - output csv expects at least one field\n");
        return INPUT_ERR_COUNT_MISMATCH;
    }

    format->csv.field_names = (char **)malloc(sizeof(char *) * (size_t)field_count);
    if (!format->csv.field_names) {
        fprintf(stderr, "Couldn't allocate memory for csv fields: %s", strerror(errno));
        return INPUT_ERR_ALLOC;
    }

    for (int i = 0; i < field_count; ++i) {
        int supported_field = 0;
        for (size_t j = 0; j < CSV_FIELD_FUNCS_COUNT; j++) {
            if (strcmp(argv[idx + i], csv_field_map[j].name) == 0) {
                supported_field = 1;
                format->csv.field_names[i] = dup_str(argv[idx + i]);
                break;
            }
        }
        if (!supported_field) {
            fprintf(stderr, "Input error - unsupported csv field %s\n", argv[idx + i]);
            return INPUT_ERR_INVALID_ARG;
        }
    }
    format->csv.field_count = field_count;
    return INPUT_ERR_NONE;
}

static InputErrorFlag parse_output_xyz(int argc, char **argv, struct LoggingState *ls)
{
    int idx = 2;

    OutputFormat *format = &(ls->out_formats[ls->out_formats_cnt - 1]);
    format->type = OUTPUT_FORMAT_XYZ;

    // Check for optional 'stripped' parameter
    if (idx < argc && strcmp(argv[idx], "stripped") == 0) {
        format->xyz.stripped = true;
        idx++;
    }

    // if prefix is present and not a schedule keyword, use it; else generate default
    if (idx < argc && strcmp(argv[idx], "interval") != 0 && strcmp(argv[idx], "list") != 0) {
        snprintf(format->xyz.prefix, sizeof(format->xyz.prefix), "%s", argv[idx]);
        idx++;
    } else {
        // generate default prefix: [time].xyz
        time_t now = time(NULL);
        snprintf(format->xyz.prefix, sizeof(format->xyz.prefix), "%ld.xyz", (long)now);
    }

    InputErrorFlag err = parse_output_schedule(argc, argv, &idx, &format->xyz.schedule, 0);
    if (err != INPUT_ERR_NONE) {
        return err;
    }
    // error if extra parameters are present (e.g., 'fields' or any other)
    if (idx < argc) {
        fprintf(stderr, "Input error - additional parameters after %s are not recognized\n",
                argv[idx]);
        return INPUT_ERR_INVALID_ARG;
    }
    // TODO: initialize suffix?
    return INPUT_ERR_NONE;
}

static InputErrorFlag parse_steps_csv(int argc, char **argv, struct LoggingState *ls)
{
    int idx = 2;

    OutputFormat *format = &(ls->out_formats[ls->out_formats_cnt - 1]);
    format->type = OUTPUT_FORMAT_STEPS_CSV;

    // if filename is present, use it; else generate default
    if (idx < argc) {
        snprintf(format->steps.filename, sizeof(format->steps.filename), "%s", argv[idx]);
        idx++;
    } else {
        // generate default filename: [time in seconds]_iter.csv
        time_t now = time(NULL);
        snprintf(format->steps.filename, sizeof(format->steps.filename), "%ld_steps.csv",
                 (long)now);
    }

    // Check for optional 'coord' parameter
    if (idx < argc && strcmp(argv[idx], "coord") == 0) {
        format->steps.with_coordination = true;
        idx++;
    }

    // error if extra parameters are present (e.g., anything else)
    if (idx < argc) {
        fprintf(stderr, "Input error - additional parameters after %s are not recognized\n",
                argv[idx]);
        return INPUT_ERR_INVALID_ARG;
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_output(int argc, char **argv, int line, ParseContext *p_ctx,
                                 struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)inputs;

    if (argc < 2) {
        fprintf(stderr, "Input error - output command expects at least 1 argument\n");
        return INPUT_ERR_COUNT_MISMATCH;
    }

    ls->out_formats =
        realloc(ls->out_formats, (size_t)(ls->out_formats_cnt + 1) * sizeof(OutputFormat));
    if (ls->out_formats == NULL) {
        fprintf(stderr, "Couldn't allocate memory for output formats: %s", strerror(errno));
        return INPUT_ERR_ALLOC;
    }

    memset(&ls->out_formats[ls->out_formats_cnt], 0, sizeof(OutputFormat));
    ls->out_formats_cnt++;

    if (strcmp(argv[1], "csv") == 0) {
        return parse_output_csv(argc, argv, ls);
    } else if (strcmp(argv[1], "xyz") == 0) {
        return parse_output_xyz(argc, argv, ls);
    } else if (strcmp(argv[1], "steps") == 0) {
        return parse_steps_csv(argc, argv, ls);
    } else {
        fprintf(stderr, "Input error - output expects 'csv' or 'xyz'\n");
        return INPUT_ERR_INVALID_ARG;
    }
}

static InputErrorFlag cmd_geometry(int argc, char **argv, int line, ParseContext *p_ctx,
                                   struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc < 3) {
        fprintf(stderr, "Input error - not enough arguments for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }
    if (strcmp(argv[1], "sheet") == 0) {
        if (argc != 3 || parse_int(argv[2], &inputs->geometry_param)) {
            fprintf(stderr, "Input error - invalid sheet thickness %s\n", argc > 2 ? argv[2] : "");
            return INPUT_ERR_INVALID_ARG;
        }
        inputs->geometry = GEOMETRY_FLAT_SHEET;
    } else if (strcmp(argv[1], "cluster") == 0) {
        if (argc != 3 || parse_int(argv[2], &inputs->geometry_param)) {
            fprintf(stderr, "Input error - invalid cluster radius %s\n", argc > 2 ? argv[2] : "");
            return INPUT_ERR_INVALID_ARG;
        }
        inputs->geometry = GEOMETRY_CLUSTER;
    } else if (strcmp(argv[1], "file") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Input error - geometry file expects a filename\n");
            return INPUT_ERR_COUNT_MISMATCH;
        }
        inputs->geometry = GEOMETRY_FROM_FILE;
        snprintf(inputs->atoms_filename, sizeof(inputs->atoms_filename), "%s", argv[2]);
    } else {
        fprintf(stderr, "Input error - unrecognized geometry type %s\n", argv[1]);
        return INPUT_ERR_INVALID_ARG;
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_run(int argc, char **argv, int line, ParseContext *p_ctx,
                              struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc != 3) {
        fprintf(stderr, "Input error - expected 2 arguments for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }
    if (strcmp(argv[1], "time") == 0) {
        inputs->sim_end_type = SIM_END_BY_STIME;
        if (parse_double(argv[2], &inputs->run_stime)) {
            fprintf(stderr, "Input error - invalid run time %s\n", argv[2]);
            return INPUT_ERR_INVALID_ARG;
        }
    } else if (strcmp(argv[1], "iteration") == 0) {
        inputs->sim_end_type = SIM_END_BY_ITERATIONS;
        if (parse_ulong(argv[2], &inputs->final_iteration)) {
            fprintf(stderr, "Input error - invalid run iteration %s\n", argv[2]);
            return INPUT_ERR_INVALID_ARG;
        }
    } else {
        fprintf(stderr, "Input error - unknown run mode %s\n", argv[1]);
        return INPUT_ERR_INVALID_ARG;
    }
    return INPUT_ERR_NONE;
}

static InputErrorFlag cmd_flavor(int argc, char **argv, int line, ParseContext *p_ctx,
                                 struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)ls;
    if (argc != 2) {
        fprintf(stderr, "Input error - expected 1 argument for command %s\n", argv[0]);
        return INPUT_ERR_COUNT_MISMATCH;
    }
    if (strcmp(argv[1], "KMC") == 0) {
        inputs->flavor = FLAVOR_KMC;
    } else if (strcmp(argv[1], "MC") == 0) {
        inputs->flavor = FLAVOR_MC;
    } else {
        fprintf(stderr, "Input error - unknown flavor %s\n", argv[1]);
        return INPUT_ERR_INVALID_ARG;
    }
    return INPUT_ERR_NONE;
}

/* ================= Help ================= */

void print_help(const char *cmd)
{
    for (const Command *c = commands; c->name != NULL; c++) {
        if (!cmd || strcmp(cmd, c->name) == 0) {
            printf("%s\n  %s\n\n", c->usage, c->description);
        }
    }
}

static InputErrorFlag cmd_help(int argc, char **argv, int line, ParseContext *p_ctx,
                               struct SimulationConfig *inputs, struct LoggingState *ls)
{
    (void)line;
    (void)p_ctx;
    (void)inputs;
    (void)ls;
    if (argc == 1)
        print_help(NULL);
    else
        print_help(argv[1]);
    return INPUT_ERR_OTHER;
}

/* ================= Command table ================= */
// TODO: need to add a way to associate lattice dimensions with cartesian dimensions
// requriements: 1 = required, 0 = optional, -1 = forbidden
const Command commands[] = {
    {"systemsize", cmd_systemsize, "systemsize NX NY NZ",
     "Simulation box size (x,y,z) in cartesian units.", CMDCAT_GEOMETRY, CMDREQ_REQUIRED,
     CMDREQ_REQUIRED},

    {"temp", cmd_temp, "temp T", "Simulation temperature in K.", CMDCAT_THERMODYNAMICS,
     CMDREQ_REQUIRED, CMDREQ_REQUIRED},

    {"seed", cmd_seed, "seed random|default|N",
     "Selection of the random seed value. Default: default", CMDCAT_RUN, CMDREQ_OPTIONAL,
     CMDREQ_OPTIONAL},

    {"potential", cmd_potential, "potential U0 [dUdt Umax]",
     "Constant or swept electric potential in eV and eV/s. Default: 0", CMDCAT_THERMODYNAMICS,
     CMDREQ_OPTIONAL, CMDREQ_FORBIDDEN},

    {"struct", cmd_struct, "struct FCC|BCC|SC", "Crystal structure type.", CMDCAT_GEOMETRY,
     CMDREQ_REQUIRED, CMDREQ_REQUIRED},

    {"output", cmd_output,
     "output [steps [filename] | csv [filename] | xyz [file prefix]] [interval|list] "
     "[time|iteration] [interval "
     "step|list] {fields [list of fields]}",
     "Data output scheme. CSV formats support {fields} fields. XYZ filename is used as a prefix, "
     "and new xyz files are output following the schedule",
     CMDCAT_OUTPUT, CMDREQ_OPTIONAL, CMDREQ_OPTIONAL},

    {"geometry", cmd_geometry, "geometry (sheet N|cluster R|file path)",
     "Initial geometry configuration. `sheet` parameter is thickness in the third lattice "
     "dimension, `cluster` parameter is the cluster radius in lattice units, and `file` parameter "
     "is a path to an input file with atomic coordinates.",
     CMDCAT_GEOMETRY, CMDREQ_REQUIRED, CMDREQ_REQUIRED},

    {"atomtype", cmd_atomtype, "atomtype A B [C ...]", "Define atom types and their order.",
     CMDCAT_GEOMETRY, CMDREQ_REQUIRED, CMDREQ_REQUIRED},

    {"composition", cmd_composition, "composition xA xB [xC ...]",
     "Atomic composition fractions; order follows atomtype.", CMDCAT_GEOMETRY, CMDREQ_REQUIRED,
     CMDREQ_REQUIRED},

    {"dissolution", cmd_dissolution, "dissolution true|false ...",
     "Dissolution flags per atom type; order follows atomtype.", CMDCAT_GEOMETRY, CMDREQ_REQUIRED,
     CMDREQ_FORBIDDEN},

    {"nnlevels", cmd_nnlevels, "nnlevels N", "Number of nearest-neighbor shells.",
     CMDCAT_THERMODYNAMICS, CMDREQ_REQUIRED, CMDREQ_REQUIRED},

    {"nne", cmd_nne, "1nne eAA eAB ...",
     "Nearest-neighbor energies for shell n (flattened upper-triangle). For a three-component "
     "system: AA AB AC BB BC CC",
     CMDCAT_THERMODYNAMICS, CMDREQ_REQUIRED, CMDREQ_REQUIRED},

    {"run", cmd_run, "run time|iteration value", "Simulation end condition, time in seconds.",
     CMDCAT_RUN, CMDREQ_REQUIRED, CMDREQ_REQUIRED},

    {"flavor", cmd_flavor, "flavor KMC|MC", "Simulation algorithm flavor.", CMDCAT_RUN,
     CMDREQ_REQUIRED, CMDREQ_REQUIRED},

    {"help", cmd_help, "help [command]", "Show documentation for commands.", CMDCAT_OUTPUT,
     CMDREQ_OPTIONAL, CMDREQ_OPTIONAL},

    {NULL, NULL, NULL, NULL, CMDCAT_UNCAT, CMDREQ_OPTIONAL, CMDREQ_OPTIONAL}};

void initialize_requirements(ParseContext *p_ctx)
{
    size_t commands_count = sizeof(commands) / sizeof(commands[0]) - 1;
    p_ctx->cmd_present = calloc(commands_count, sizeof(int));
    if (p_ctx->cmd_present == NULL) {
        fprintf(stderr, "Couldn't allocate memory for command requirements in ParseContext: %s",
                strerror(errno));
        call_exit(INPUT_ERR_ALLOC);
    }
}

void mark_requirement(ParseContext *p_ctx, const Command *c)
{
    ptrdiff_t idx = c - commands;
    p_ctx->cmd_present[idx] = 1;
}

/* ================= Parser ================= */

void clean_ctx(ParseContext *p_ctx)
{
    for (int i = 0; i < p_ctx->nne_cmd_count; i++) {
        free_if_exists((void **)&p_ctx->nne_cmds[i].values);
    }
    free_if_exists((void **)&p_ctx->nne_cmds);
    free_if_exists((void **)&p_ctx->composition_raw);
    free_if_exists((void **)&p_ctx->dissolution_raw);
    free_if_exists((void **)&p_ctx->cmd_present);
}

void parse_input_file(FILE *fp, ParseContext *p_ctx, struct SimulationConfig *inputs,
                      struct LoggingState *ls)
{
    char line[MAX_LINE];
    char *argv[MAX_TOKENS];
    int lineno = 0;

    initialize_requirements(p_ctx);
    while (fgets(line, sizeof(line), fp)) {
        lineno++;

        // ignore comments
        if (line[0] == '#')
            continue;

        int argc = tokenize_line(line, argv, MAX_TOKENS);
        if (argc == 0)
            continue;

        int matched = 0;
        for (const Command *c = commands; c->name != NULL; c++) {
            if (strcmp(argv[0], c->name) == 0 ||
                (strcmp(c->name, "nne") == 0 && isdigit(argv[0][0]))) {
                mark_requirement(p_ctx, c);
                InputErrorFlag err = c->func(argc, argv, lineno, p_ctx, inputs, ls);
                if (err != INPUT_ERR_NONE) {
                    // if handler fails
                    fprintf(stderr, "Error in command '%s' at line %d\nUsage: %s\n", argv[0],
                            lineno, c->usage);
                    call_exit(err);
                }
                matched = 1;
                break;
            }
        }

        if (!matched) {
            fprintf(stderr, "Unknown command '%s' at line %d\n", argv[0], lineno);
            call_exit(INPUT_ERR_UNKNOWN_CMD);
        }
    }
}

/* ================= Finalization ================= */

void finalize_atom_dependent(const ParseContext *p_ctx, struct SimulationConfig *inputs)
{
    if (inputs->atom_names_cnt <= 0) {
        fprintf(stderr, "atomtype must be specified\n");
        clean_and_error(INPUT_ERR_MISSING_DEPENDENCY);
    }

    if (p_ctx->composition_raw) {
        if (p_ctx->composition_count != inputs->atom_names_cnt) {
            fprintf(stderr, "composition count mismatch at line %d\n", p_ctx->composition_lineno);
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_COUNT_MISMATCH);
        }
        double total = 0.0;
        for (int i = 0; i < p_ctx->composition_count; i++) {
            total += p_ctx->composition_raw[i];
        }
        if (fabs(total - 1.0) > 1e-10) {
            fprintf(stderr, "composition values must sum to 1.0, got %.17g\n", total);
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_MISSING_DEPENDENCY);
        }
        inputs->substrate_composition =
            (double *)malloc((size_t)inputs->atom_names_cnt * sizeof(double));
        if (inputs->substrate_composition == NULL) {
            fprintf(stderr, "Couldn't allocate memory for substrate composition: %s",
                    strerror(errno));
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_ALLOC);
        }
        memcpy(inputs->substrate_composition, p_ctx->composition_raw,
               (size_t)inputs->atom_names_cnt * sizeof(double));
    }

    // initialize solubility flags, default to false if not specified
    inputs->is_soluble = (bool *)malloc((size_t)inputs->atom_names_cnt * sizeof(bool));
    if (inputs->is_soluble == NULL) {
        fprintf(stderr, "Couldn't allocate memory for dissolution flags, is_soluble: %s",
                strerror(errno));
        clean_ctx((ParseContext *)p_ctx);
        clean_and_error(INPUT_ERR_ALLOC);
    }

    inputs->dissolution = 0;
    if (p_ctx->dissolution_raw) {
        if (p_ctx->dissolution_count != inputs->atom_names_cnt) {
            fprintf(stderr, "dissolution count mismatch at line %d\n", p_ctx->dissolution_lineno);
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_COUNT_MISMATCH);
        }
        memcpy(inputs->is_soluble, p_ctx->dissolution_raw,
               (size_t)inputs->atom_names_cnt * sizeof(bool));

        for (int i = 0; i < p_ctx->dissolution_count; i++) {
            if (p_ctx->dissolution_raw[i]) {
                inputs->dissolution = 1;
                break;
            }
        }
    } else {
        for (int i = 0; i < inputs->atom_names_cnt; i++) {
            inputs->is_soluble[i] = false;
        }
    }
}

void finalize_nne(const ParseContext *p_ctx, struct SimulationConfig *inputs)
{
    int ntypes = inputs->atom_names_cnt;
    inputs->num_bond_types = get_num_bond_types(inputs->num_elements);

    if (inputs->num_nn_levels <= 0) {
        fprintf(stderr, "nnlevels must be specified\n");
        clean_ctx((ParseContext *)p_ctx);
        clean_and_error(INPUT_ERR_MISSING_DEPENDENCY);
    }

    inputs->num_nn_types = inputs->num_nn_levels * inputs->num_bond_types;
    inputs->nn_energy = (double *)calloc((size_t)inputs->num_nn_types, sizeof(double));

    int *defined = calloc((size_t)inputs->num_nn_levels, sizeof(int));
    int expected = inputs->num_bond_types;

    for (int c = 0; c < p_ctx->nne_cmd_count; c++) {
        int lvl = p_ctx->nne_cmds[c].level;
        if (p_ctx->nne_cmds[c].nvalues != expected) {
            fprintf(stderr, "nne at line %d expects %d values, got %d\n", p_ctx->nne_cmds[c].line,
                    expected, p_ctx->nne_cmds[c].nvalues);
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_COUNT_MISMATCH);
        }
        if (lvl > inputs->num_nn_levels) {
            fprintf(stderr, "nne level %d exceeds nnlevels (line %d)\n", lvl,
                    p_ctx->nne_cmds[c].line);
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_INVALID_ARG);
        }
        if (defined[lvl - 1]) {
            fprintf(stderr, "duplicate nne for level %d (line %d)\n", lvl, p_ctx->nne_cmds[c].line);
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_DUPLICATE_CMD);
        }

        if (p_ctx->nne_cmds[c].values == NULL) {
            fprintf(stderr, "nne command at line %d has no values\n", p_ctx->nne_cmds[c].line);
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_COUNT_MISMATCH);
        }

        // Fill upper triangle from flattened input
        int bond_index = nn_bondidx_2_envidx(lvl - 1, 0, inputs->num_bond_types);
        double *vals = p_ctx->nne_cmds[c].values;
        int idx = 0;
        for (int i = 0; i < ntypes; i++) {
            for (int j = i; j < ntypes; j++) {
                if (idx >= expected)
                    break;
                inputs->nn_energy[bond_index] = vals[idx];
                idx++;
                bond_index++;
            }
        }
        defined[lvl - 1] = 1;
    }

    for (int i = 0; i < inputs->num_nn_levels; i++) {
        if (!defined[i]) {
            fprintf(stderr, "missing nne definition for level %d\n", i + 1);
            clean_ctx((ParseContext *)p_ctx);
            ;
            clean_and_error(INPUT_ERR_MISSING_CMD);
        }
    }

    free(defined);
}

// check that all fields specified in csv output are supported by the specified flavor
// populate field_funcs arrays based on names of fields
void finalize_csv_fields(struct LoggingState *ls, unsigned int flavor)
{

    for (int o = 0; o < ls->out_formats_cnt; o++) {
        OutputFormat *format = &ls->out_formats[o];
        if (format->type != OUTPUT_FORMAT_CSV) {
            continue;
        }

        format->csv.field_funcs =
            (CsvFieldFuncPtr *)malloc(sizeof(CsvFieldFuncPtr *) * (size_t)format->csv.field_count);
        if (!format->csv.field_funcs) {
            fprintf(stderr, "Couldn't allocate memory for csv fields: %s", strerror(errno));
            clean_and_error(INPUT_ERR_ALLOC);
        }

        char *flavor_name = flavor == 1 ? "KMC" : (flavor == 2 ? "MC" : "Undefined");

        // iterate over user-specified csv fields
        for (int i = 0; i < format->csv.field_count; ++i) {
            const char *field_name = format->csv.field_names[i];

            // check field name against supported fields
            for (size_t j = 0; j < CSV_FIELD_FUNCS_COUNT; j++) {
                if (strcmp(field_name, csv_field_map[j].name) == 0) {
                    // if generally supported, check if supported for this flavor
                    if ((csv_field_map[j].flavor != FLAVOR_UNDEFINED) &&
                        (csv_field_map[j].flavor != flavor)) {
                        fprintf(stderr, "Input error - field '%s' is not supported for flavor %s\n",
                                field_name, flavor_name);
                        clean_and_error(INPUT_ERR_INVALID_ARG);
                    }
                    format->csv.field_funcs[i] = csv_field_map[j].get_value;
                    break;
                }
            }
        }
    }
}

void check_required_inputs(const ParseContext *p_ctx, unsigned int flavor)
{
    int failed_requirements = 0;
    size_t commands_count = sizeof(commands) / sizeof(commands[0]) - 1;
    for (size_t i = 0; i < commands_count; i++) {
        // =1 if command is required and not present
        // =0 if command is required and present, or not required
        // =-1 if command is present and forbidden/ignored
        int present = p_ctx->cmd_present[i];
        int req;
        if (flavor == FLAVOR_MC) {
            req = commands[i].required_mc;
        } else if (flavor == FLAVOR_KMC) {
            req = commands[i].required_kmc;
        } else {
            fprintf(stderr, "Unknown simulation flavor %u\n", flavor);
            clean_ctx((ParseContext *)p_ctx);
            clean_and_error(INPUT_ERR_MISSING_CMD);
        }

        // ok: required && present, optional && !present, optional && present, forbidden && !present
        if ((req == CMDREQ_REQUIRED) && !present) {
            fprintf(stderr, "Required command %s was not provided\n", commands[i].name);
            failed_requirements = 1;
        } else if ((req == CMDREQ_FORBIDDEN) && present) {
            fprintf(stderr, "Command %s is not applicable for flavor %u - ignored\n",
                    commands[i].name, flavor);
        }
    }
    if (failed_requirements) {
        clean_and_error(INPUT_ERR_MISSING_CMD);
    }
}

void validate_input_values(const struct SimulationConfig *inputs)
{
    if (inputs->system_size_x <= 0 || inputs->system_size_y <= 0 || inputs->system_size_z <= 0) {
        fprintf(stderr, "System size must be positive in all dimensions\n");
        clean_and_error(INPUT_ERR_INVALID_ARG);
    }
    if (inputs->temperature <= 0) {
        fprintf(stderr, "Temperature must be positive\n");
        clean_and_error(INPUT_ERR_INVALID_ARG);
    }
    if (inputs->num_nn_levels <= 0) {
        fprintf(stderr, "Number of nearest neighbor levels must be positive\n");
        clean_and_error(INPUT_ERR_INVALID_ARG);
    } else if (inputs->num_nn_levels > MAX_NN_LEVELS) {
        fprintf(stderr, "Number of nearest neighbor levels must not exceed %d\n", MAX_NN_LEVELS);
        clean_and_error(INPUT_ERR_INVALID_ARG);
    }
}

void finalize_config(const ParseContext *p_ctx, struct SimulationConfig *inputs,
                     struct LoggingState *ls)
{
    finalize_atom_dependent(p_ctx, inputs);
    finalize_nne(p_ctx, inputs);
    finalize_csv_fields(ls, inputs->flavor);
    check_required_inputs(p_ctx, inputs->flavor);
    validate_input_values(inputs);
}

// finalize_nne using multidimensional array for nn_energy

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
