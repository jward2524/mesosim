#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "Input.h"
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
#endif

static void usage(int help_type)
{
    if (help_type == 2) {
        printf("\nInput file commands:\n\n");
        print_help(NULL);
    } else {
        printf("Usage: mesosim [OPTIONS] [FILE]\n"
               "Execute the mesosim KMC simulation program, with FILE as the input file\n"
               "\n"
               "Options:\n"
               "  -h, --help\tDisplay this message\n"
               "  -i, --help-input\tDisplay input file help\n"
               "  -v, --verbose\tPrint additional information\n");
    }
    call_exit(EXIT_SUCCESS);
}

static void parse_arguments(int argc, char *argv[], char **pfilename, int *verbose_flag,
                            unsigned long *verbose_interval, int *help_type)
{
    if (argc <= 1) {
        pfilename = NULL;
        return;
    }

    *pfilename = NULL;
    *verbose_flag = TEST;

    // TODO: add a flag for starting from a checkpoint file
    // argv[0] is executable / program name
    for (int i = 1; i < argc; i++) {
        char *arg = argv[i];
        if ((strcmp(arg, "-h") == 0) || (strcmp(arg, "--help") == 0)) {
            *pfilename = NULL;
            *help_type = 1;
            return;
        } else if ((strcmp(arg, "-i") == 0) || (strcmp(arg, "--help-input") == 0)) {
            *pfilename = NULL;
            *help_type = 2;
            return;
        } else if ((strcmp(arg, "-v") == 0) || (strcmp(arg, "--verbose") == 0)) {
            *verbose_flag = 1;
            *verbose_interval = strtoul(argv[i + 1], NULL, 10);
            if (*verbose_interval == 0) {
                *verbose_interval = 1000;
            } else {
                i++;
            }
        } else if (arg[0] == '-') {
            // if it starts with a tack, assume it is an option that can't be parsed
            fprintf(stderr, "Unrecognized option: %s\n\n", arg);
            *pfilename = NULL;
            return;
        } else {
            // if no tack, assume it is the filename
            if (*pfilename) {
                fprintf(stderr,
                        "Only one input file is supported: %s is the second unsupported argument\n",
                        arg);
            } else {
                *pfilename = arg;
            }
        }
    }
}

static void simulation_setup(struct SimulationState *ss, struct SimulationEnv *se,
                             struct LoggingState *ls);
static void perform_simulation(struct SimulationState *ss, struct SimulationEnv *se,
                               struct LoggingState *ls);
static void simulation_cleanup(struct SimulationState *ss, struct SimulationEnv *se,
                               struct LoggingState *ls);

int main(int argc, char *argv[])
{
    clock_t presim_start = clock();
    char *input_filename = NULL;
    int verbose_flag = 0;
    unsigned long verbose_interval = 0;
    int help_type = 0;
    parse_arguments(argc, argv, &input_filename, &verbose_flag, &verbose_interval, &help_type);
    if (!input_filename) {
        usage(help_type);
    }

    // start the time
    time_t starttime = 0;
    time(&starttime);

    struct SimulationState *sim_state;
    struct SimulationEnv *sim_env;
    struct LoggingState *log_state;
    initialize_states(&sim_state, &sim_env, &log_state);

    log_state->sim_log = stdout;
    // fresh start here
    log_state->verbose = verbose_flag;
    log_state->verbose_interval = verbose_interval;

    safe_log(log_state->sim_log, "Start time: %s\n", ctime(&starttime));
    safe_log(log_state->sim_log, "Attempting to read in file %s\n", input_filename);

    // pre-process the file information and fill in the gaps with defaults
    struct SimulationConfig inputs = {0};
    simulation_parameters_from_file(input_filename, &inputs, log_state);
    printf("Read file successfully\n");

    open_log_files(log_state, inputs.flavor);

    initialize_simulation(&inputs, sim_state, sim_env, log_state);
    // end fresh start

    double presim_time = (double)(clock() - presim_start) / CLOCKS_PER_SEC;
    safe_log(log_state->sim_log, "Pre-simulation setup time: %lg seconds\n", presim_time);

    clock_t sim_start = clock();

    simulation_setup(sim_state, sim_env, log_state);
    perform_simulation(sim_state, sim_env, log_state);
    simulation_cleanup(sim_state, sim_env, log_state);

    double sim_time = (double)(clock() - sim_start) / CLOCKS_PER_SEC;
    safe_log(log_state->sim_log, "Simulation time: %lg seconds\n", sim_time);
    safe_log(log_state->sim_log, "Average iteration time: %lg seconds\n",
             sim_time / (double)sim_state->iter);

    clean_and_error(EXIT_SUCCESS);

    return 0;
}

static void simulation_setup(struct SimulationState *ss, struct SimulationEnv *se,
                             struct LoggingState *ls)
{
    if (ls->verbose) {
        printf("Beginning simulation\n");
    }

    // refresh transitions and compute transition array before starting simulation loop
    for (int i = 0; i < ss->atom_cnt; ++i) {
        refresh_transitions(i, ss, se);
    }
    if (se->flavor == FLAVOR_KMC) {
        compute_transition_array(ss, se);
    }

    // log initial state at beginning of simulation
    // if starting from checkpoint, don't log initial state bc already logged
    if (ss->iter == 0) {
        write_logs(NULL, ss, se, ls);
    }

    return;
}

static void perform_simulation(struct SimulationState *ss, struct SimulationEnv *se,
                               struct LoggingState *ls)
{
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
