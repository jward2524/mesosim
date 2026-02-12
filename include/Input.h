#ifndef KMC_INPUT_H
#define KMC_INPUT_H

#include "State.h"
#include <stdio.h>

/* ================= Simulation configuration ================= */

typedef struct {
    int natomtypes;
    char **atomtypes;

    double *composition;  // size natomtypes
    int *dissolution;     // size natomtypes

    int nnlevels;
    double ***nne;        // nne[shell][i][j], symmetric

} SimConfig;

/* ================= Deferred parse context ================= */

struct NNECmd {
    int level;
    double *values;   // flattened upper-triangle for arbitrary ntypes
    int nvalues;
    int line;
};

typedef struct {
    struct NNECmd *nne_cmds;
    int nne_cmd_count;

    double *composition_raw;
    int composition_count;
    int composition_line;

    bool *dissolution_raw;
    int dissolution_count;
    int dissolution_line;
} ParseContext;

/* ================= Command Table =============== */

enum CmdCategory {
    CMDCAT_UNCAT = -1,
    CMDCAT_GEOMETRY = 0,
    CMDCAT_THERMODYNAMICS,
    CMDCAT_OUTPUT,
    CMDCAT_RUN,
};

typedef int (*CmdFunc)(int argc, char **argv, int line, ParseContext *ctx,
                       struct SimulationState *ss, struct SimulationEnv *se,
                       struct LoggingState *ls);

typedef struct {
    const char *name;
    CmdFunc func;
    const char *usage;
    const char *description;
    const enum CmdCategory category;
    const int required;
} Command;

extern Command commands[];

/* ================= Parsing API ================= */

void parse_input_file(FILE *fp, ParseContext *ctx, struct SimulationState *ss,
                      struct SimulationEnv *se, struct LoggingState *ls);
void finalize_atom_dependent(const ParseContext *ctx, struct SimulationEnv *se);
void finalize_nne(const ParseContext *ctx, struct SimulationEnv *se);
void finalize_config(const ParseContext *ctx, struct SimulationState *ss,
                     struct SimulationEnv *se, struct LoggingState *ls);
void print_help(const char *cmd);

#endif
