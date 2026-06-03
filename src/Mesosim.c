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

static struct SimulationState *sim_state;
static struct SimulationEnv *sim_env;
static struct LoggingState *log_state;

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

    initialize_states(&sim_state, &sim_env, &log_state);

    log_state->verbose = verbose_flag;
    log_state->verbose_interval = verbose_interval;
    log_state->sim_log = stdout;

    safe_log(log_state->sim_log, "Start time: %s\n", ctime(&starttime));
    safe_log(log_state->sim_log, "Attempting to read in file %s\n", input_filename);

    // pre-process the file information and fill in the gaps with defaults
    struct SimulationConfig inputs = {0};
    simulation_parameters_from_file(input_filename, &inputs, log_state);
    printf("Read file successfully\n");

    open_log_files(log_state, inputs.flavor);

    // TODO: take user inputs as argument
    initialize_simulation(&inputs, sim_state, sim_env, log_state);

    double presim_time = (double)(clock() - presim_start) / CLOCKS_PER_SEC;
    safe_log(log_state->sim_log, "Pre-simulation setup time: %lg seconds\n", presim_time);

    if (log_state->verbose)
        printf("Beginning simulation\n");

    clock_t sim_start = clock();

    // perform simulations
    unsigned long sim_error;
    switch (sim_env->flavor) {
    case FLAVOR_KMC:
        sim_error = perform_simulation(sim_state, sim_env, log_state);
        break;
    case FLAVOR_MC:
        sim_error = perform_metropolis_mc(sim_state, sim_env, log_state);
        break;
    default:
        fprintf(stderr, "Flavor %d not recognized, no simulation performed\n", sim_env->flavor);
        sim_error = 1;
        break;
    }

    if (sim_error != 0) {
        fprintf(stderr, "ERROR! Something went wrong in the simulation\n");
        clean_and_error(EXIT_FAILURE);
    }
    double sim_time = (double)(clock() - sim_start) / CLOCKS_PER_SEC;
    safe_log(log_state->sim_log, "Simulation time: %lg seconds\n", sim_time);
    safe_log(log_state->sim_log, "Average iteration time: %lg seconds\n",
             sim_time / (double)sim_state->iter);

    clean_and_error(EXIT_SUCCESS);

    return 0;
}
