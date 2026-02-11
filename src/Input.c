#include "Input.h"
#include "Utils.h"
#include "ErrorM.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE   1024
#define MAX_TOKENS 256

// if function returns 0, function failed
// if it returns 1, handler succeeeded

/* ================= Utilities ================= */
// TODO: replace with XYZParser.c:tokenize_line
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
    long v = strtol(s, &e, 10);
    if (e == NULL) {
        return 1;
    }
    *out = (int)v;
    return 0;
}

static int parse_double(const char *s, double *out) {
    char *e;
    double v = strtod(s, &e);
    if (e == NULL) {
        return 1;
    }
    *out = v;
    return 0;
}

/* ================= Command handlers ================= */
typedef int (*CmdFunc)(int argc, char **argv, int line, ParseContext *ctx,
                       struct SimulationState *ss, struct SimulationEnv *se,
                       struct LoggingState *ls);

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
    if (argc < 2) {
        fprintf(stderr, "Input error - no arguments provided for command %s\n", argv[0]);
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
        values[i] = strtod(argv[i+1], NULL);
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

/* ================= Command table ================= */
static int cmd_help(int argc, char **argv, int line, ParseContext *ctx,
                    struct SimulationState *ss, struct SimulationEnv *se,
                    struct LoggingState *ls);

typedef struct {
    const char *name;
    CmdFunc func;
    const char *usage;
    const char *description;
} Command;


static Command commands[] = {
    { "atomtype",    cmd_atomtype,    "atomtype A B [C ...]",  "Define atom types and their order." },
    { "composition", cmd_composition, "composition xA xB [xC ...]", "Atomic composition fractions; order follows atomtype." },
    { "dissolution", cmd_dissolution, "dissolution true|false ...", "Dissolution flags per atom type; order follows atomtype." },
    { "nnlevels",    cmd_nnlevels,    "nnlevels N", "Number of nearest-neighbor shells." },
    { "nne",         cmd_nne,         "1nne eAA eAB ...", "Nearest-neighbor energies for shell n (flattened upper-triangle)." },
    { "help",        cmd_help,        "help [command]", "Show documentation for commands." },
    { NULL, NULL, NULL, NULL }
};

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


// TODO: move to Input.c - currently will crash on build
/* ================= Finalization ================= */

static void finalize_atom_dependent(const ParseContext *ctx, struct SimulationEnv *se)
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
    }
}

static void finalize_nne(const ParseContext *ctx, struct SimulationEnv *se)
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

    finalize_config(ctx, ss, se, ls);
}
