#ifndef KMC_INPUT_H
#define KMC_INPUT_H

#include "State.h"
#include <stdio.h>

/* ================= Input Error Flags ================= */
/**
 * Error flags for input parsing and validation.
 * Used to distinguish error types in call_exit/clean_and_error.
 */
typedef enum {
    INPUT_ERR_NONE = 0,           // No error
    INPUT_ERR_INVALID_ARG,        // Invalid argument or value
    INPUT_ERR_MISSING_CMD,        // Missing required command
    INPUT_ERR_DUPLICATE_CMD,      // Duplicate command
    INPUT_ERR_COUNT_MISMATCH,     // Argument count mismatch
    INPUT_ERR_MISSING_DEPENDENCY, // Command dependency was not met
    INPUT_ERR_UNKNOWN_CMD,        // Unknown command
    INPUT_ERR_ALLOC,              // Allocation failure
    INPUT_ERR_OTHER               // Other/unclassified error
} InputErrorFlag;

/* ================= Simulation configuration ================= */

typedef struct {
    int natomtypes;
    char **atomtypes;

    double *composition; // size natomtypes
    int *dissolution;    // size natomtypes

    int nnlevels;
    double ***nne; // nne[shell][i][j], symmetric

} SimConfig;

/* ================= Deferred parse context ================= */

struct NNECmd {
    int level;
    double *values; // flattened upper-triangle for arbitrary ntypes
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

    // array that holds flags for whether required commands were provided
    int *cmd_present;
} ParseContext;

/* ================= Command Table =============== */

enum CmdCategory {
    CMDCAT_UNCAT = -1,
    CMDCAT_GEOMETRY = 0,
    CMDCAT_THERMODYNAMICS,
    CMDCAT_OUTPUT,
    CMDCAT_RUN,
};

enum CmdRequired {
    CMDREQ_REQUIRED,
    CMDREQ_OPTIONAL,
    CMDREQ_FORBIDDEN,
};

typedef InputErrorFlag (*CmdFunc)(int argc, char **argv, int line, ParseContext *ctx,
                                  struct SimulationState *ss, struct SimulationEnv *se,
                                  struct LoggingState *ls);

typedef struct {
    const char *name;
    CmdFunc func;
    const char *usage;
    const char *description;
    const enum CmdCategory category;
    const enum CmdRequired required_kmc;
    const enum CmdRequired required_mc;
} Command;

extern const Command commands[];

/* ================= Parsing API ================= */

void parse_input_file(FILE *fp, ParseContext *ctx, struct SimulationState *ss,
                      struct SimulationEnv *se, struct LoggingState *ls);
void finalize_atom_dependent(const ParseContext *ctx, struct SimulationEnv *se);
void finalize_nne(const ParseContext *ctx, struct SimulationEnv *se);
void finalize_csv_fields(struct LoggingState *ls, unsigned int flavor);
void check_required_inputs(const ParseContext *ctx, unsigned int flavor);
void finalize_config(const ParseContext *ctx, struct SimulationState *ss, struct SimulationEnv *se,
                     struct LoggingState *ls);
void print_help(const char *cmd);

#endif
