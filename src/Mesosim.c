#include "Checkpoint.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "Input.h"
#include "KMC.h"
#include "MC.h"
#include "Simulation.h"
#include "State.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// makefile defines TEST macro in debug mode
#ifndef TEST
#define TEST 0
#define MESOSIM_MAIN main
#else
#define MESOSIM_MAIN main_mesosim
#endif

typedef enum { START_INPUT = 1, START_CHECKPOINT } StartType;
typedef enum { HELP_NONE = 0, HELP_GENERAL = 1, HELP_INPUT = 2 } HelpType;

typedef struct {
    StartType start_type;
    char *filename;
    int verbose_flag;
    unsigned long verbose_interval;
    HelpType help_type;
} CmdArgs;

static void usage(int help_type);
static void parse_arguments(int argc, char *argv[], CmdArgs *args);
static void initialize_simulation(const CmdArgs cmd_args, struct SimulationState *ss,
                                  struct SimulationEnv *se, struct LoggingState *ls);
static void initialize_kinetics(struct SimulationState *ss, struct SimulationEnv *se);
static void perform_simulation(struct SimulationState *ss, struct SimulationEnv *se,
                               struct LoggingState *ls);
static void simulation_cleanup(struct SimulationState *ss, struct SimulationEnv *se,
                               struct LoggingState *ls);

int MESOSIM_MAIN(int argc, char *argv[])
{
    // simulation fails fast, so main won't return a value unless simulation successfully completes

    clock_t presim_start = clock();
    CmdArgs cmd_args = {0};
    parse_arguments(argc, argv, &cmd_args);
    if ((cmd_args.help_type != HELP_NONE) || !cmd_args.filename) {
        usage(cmd_args.help_type);
        return 1;
    }

    time_t starttime = 0;
    time(&starttime);

    struct SimulationState *sim_state;
    struct SimulationEnv *sim_env;
    struct LoggingState *log_state;
    initialize_states(&sim_state, &sim_env, &log_state);

    log_state->sim_log = stdout;
    safe_log(log_state->sim_log, "Start time: %s\n", ctime(&starttime));

    initialize_simulation(cmd_args, sim_state, sim_env, log_state);

    double presim_time = (double)(clock() - presim_start) / CLOCKS_PER_SEC;
    safe_log(log_state->sim_log, "Simulation setup time: %lg seconds\n", presim_time);

    clock_t sim_start = clock();
    perform_simulation(sim_state, sim_env, log_state);
    simulation_cleanup(sim_state, sim_env, log_state);

    double sim_time = (double)(clock() - sim_start) / CLOCKS_PER_SEC;

    safe_log(log_state->sim_log, "Simulation time: %lg seconds\n", sim_time);
    safe_log(log_state->sim_log, "Average iteration time: %lg seconds\n",
             sim_time / (double)sim_state->iter);

    clean_and_error(EXIT_SUCCESS);

    return 0;
}

static void usage(int help_type)
{
    switch (help_type) {
    case HELP_NONE:
        break;
    case HELP_GENERAL:
        printf("Usage: mesosim [OPTIONS] [FILE]\n"
               "Execute the mesosim KMC/MC simulation program, with FILE as the input file\n"
               "\n"
               "Options:\n"
               "  -h, --help\tDisplay this message\n"
               "  -i, --help-input\tDisplay input file help\n"
               "  -v, --verbose\tPrint additional information\n"
               "  -c, --checkpoint [FILE]\tStart simulation from checkpoint file [FILE]\n");
        break;
    case HELP_INPUT:
        printf("\nInput file commands:\n\n");
        print_help(NULL);
        break;
    default:
        fprintf(stderr, "Invalid help type: %d\n", help_type);
        break;
    }
}

static void parse_arguments(int argc, char *argv[], CmdArgs *args)
{
    if (argc <= 1) {
        args->filename = NULL;
        return;
    }

    args->filename = NULL;
    args->verbose_flag = TEST;

    // argv[0] is executable / program name
    for (int i = 1; i < argc; i++) {
        char *arg = argv[i];
        if ((strcmp(arg, "-h") == 0) || (strcmp(arg, "--help") == 0)) {
            args->filename = NULL;
            args->help_type = HELP_GENERAL;
            return;
        } else if ((strcmp(arg, "-i") == 0) || (strcmp(arg, "--help-input") == 0)) {
            args->filename = NULL;
            args->help_type = HELP_INPUT;
            return;
        } else if ((strcmp(arg, "-v") == 0) || (strcmp(arg, "--verbose") == 0)) {
            args->verbose_flag = 1;
            args->verbose_interval = strtoul(argv[i + 1], NULL, 10);
            if (args->verbose_interval == 0) {
                args->verbose_interval = 1000;
            } else {
                i++;
            }
        } else if ((strcmp(arg, "-c") == 0) || (strcmp(arg, "--checkpoint") == 0)) {
            if (args->filename) {
                fprintf(stderr,
                        "Only one input file is supported: %s %s is the second unsupported file\n",
                        argv[i], argv[i + 1]);
                args->help_type = HELP_GENERAL;
                return;
            }
            i++;
            if (i >= argc) {
                fprintf(stderr, "Checkpoint option provided without checkpoint file\n");
                args->help_type = HELP_GENERAL;
                return;
            }
            args->filename = argv[i];
            args->start_type = START_CHECKPOINT;
        } else if (arg[0] == '-') {
            // if it starts with a tack, assume it is an option that can't be parsed
            fprintf(stderr, "Unrecognized option: %s\n\n", arg);
            args->filename = NULL;
            args->help_type = HELP_GENERAL;
            return;
        } else {
            // if no tack, assume it is the filename
            if (args->filename) {
                fprintf(stderr,
                        "Only one input file is supported: %s is the second unsupported file\n",
                        arg);
                args->help_type = HELP_GENERAL;
                return;
            }
            args->filename = arg;
            args->start_type = START_INPUT;
        }
    }
}

static void initialize_simulation(const CmdArgs cmd_args, struct SimulationState *ss,
                                  struct SimulationEnv *se, struct LoggingState *ls)
{
    if (cmd_args.start_type == START_INPUT) {
        ls->verbose = cmd_args.verbose_flag;
        ls->verbose_interval = cmd_args.verbose_interval;

        safe_log(ls->sim_log, "Reading input file %s\n", cmd_args.filename);

        // pre-process the file information and fill in the gaps with defaults
        struct SimulationConfig inputs = {0};
        simulation_parameters_from_file(cmd_args.filename, &inputs, ls);
        printf("Read file successfully\n");

        open_log_files(ls, inputs.flavor);

        initialize_simulation_from_input(&inputs, ss, se, ls);

    } else if (cmd_args.start_type == START_CHECKPOINT) {
        safe_log(ls->sim_log, "Reading checkpoint file %s\n", cmd_args.filename);

        CheckpointStatus status = read_checkpoint_file(cmd_args.filename, ss, se, ls);

        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Failed to read checkpoint file %s\n", cmd_args.filename);
            clean_and_error(EXIT_FAILURE);
            return;
        }

    } else {
        fprintf(stderr, "Invalid start type: %d\n", cmd_args.start_type);
        clean_and_error(EXIT_FAILURE);
        return;
    }

    initialize_kinetics(ss, se);

    // log initial state at beginning of simulation
    // if starting from checkpoint, don't log initial state bc already logged
    if (cmd_args.start_type == START_INPUT) {
        write_logs(NULL, ss, se, ls);
    }
}

static void initialize_kinetics(struct SimulationState *ss, struct SimulationEnv *se)
{
    // refresh transitions and compute transition array before starting simulation loop
    for (int i = 0; i < ss->atom_cnt; ++i) {
        refresh_transitions(i, ss, se);
    }

    if (se->flavor == FLAVOR_KMC) {
        compute_transition_array(ss, se);
    }

    return;
}

static void perform_simulation(struct SimulationState *ss, struct SimulationEnv *se,
                               struct LoggingState *ls)
{
    if (ls->verbose) {
        printf("Beginning simulation\n");
    }

    unsigned long sim_error;
    switch (se->flavor) {
    case FLAVOR_KMC:
        sim_error = perform_kmc(ss, se, ls);
        break;
    case FLAVOR_MC:
        sim_error = perform_metropolis_mc(ss, se, ls);
        break;
    default:
        fprintf(stderr, "Flavor %d not recognized, no simulation performed\n", se->flavor);
        sim_error = 1;
        break;
    }
}

static void simulation_cleanup(struct SimulationState *ss, struct SimulationEnv *se,
                               struct LoggingState *ls)
{
    write_logs(NULL, ss, se, ls);

    switch (se->flavor) {
    case FLAVOR_KMC:
        if ((ss->final_iteration > 0) && (ss->iter >= ss->final_iteration)) {
            safe_log(ls->sim_log, "Reached final iteration and terminated\n");
        }
        if ((ss->run_stime > 0) && (ss->elapsed_stime >= ss->run_stime)) {
            safe_log(ls->sim_log, "Reached end of simulation time and terminated\n");
        }
        break;
    case FLAVOR_MC:
        safe_log(ls->sim_log, "Reached final iteration and terminated\n");
        break;
    }
    if (ls->verbose) {
        printf("Finished simulation\n");
    }
}
