#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
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

static time_t starttime = 0;
static time_t endtime = 0;

static struct SimulationState *sim_state;
static struct SimulationEnv *sim_env;
static struct LoggingState *log_state;

static void usage(void)
{
    printf("Usage: mesosim [OPTIONS] [FILE]\n"
           "Execute the mesosim KMC simulation program, with FILE as the input file\n"
           "\n"
           "Options:\n"
           "  -h, --help\tDisplay this message\n"
           "  -v, --verbose\tPrint additional information\n");
    call_exit(EXIT_SUCCESS);
}

static void parse_arguments(int argc, char *argv[], char **pfilename, int *verbose_flag)
{
    if (argc <= 1) {
        pfilename = NULL;
        return;
    }

    *pfilename = NULL;
    *verbose_flag = TEST;

    // argv[0] is executable / program name
    for (int i = 1; i < argc; i++) {
        char *arg = argv[i];
        if ((strcmp(arg, "-h") == 0) || (strcmp(arg, "--help") == 0)) {
            *pfilename = NULL;
            return;
        } else if ((strcmp(arg, "-v") == 0) || (strcmp(arg, "--verbose") == 0)) {
            *verbose_flag = 1;
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

    char *input_filename = NULL;
    int verbose_flag = 0;
    parse_arguments(argc, argv, &input_filename, &verbose_flag);
    if (!input_filename) {
        usage();
    }

    // start the time
    time(&starttime);

    initialize_states(&sim_state, &sim_env, &log_state);

    log_state->verbose = verbose_flag;

    // write to a temporary file until a logfile is identified
    //  not supported for msvcrt.dll [msys's mingw64]
    FILE *temp_log = tmpfile();
    if (temp_log == NULL) {
        perror("Failed to make temp file");
        clean_and_error(errno);
    }
    if (log_state->verbose)
        printf("Temporary log created\n");
    fputs("MESOSIM 2024\n", temp_log);
    safe_log(temp_log, "Start time: %s\n", ctime(&starttime));
    safe_log(temp_log, "Attempting to read in file %s\n", input_filename);

    // pre-process the file information and fill in the gaps with defaults
    bool res = simulation_parameters_from_file(input_filename, sim_state, sim_env, log_state,
                                               temp_log, starttime);

    int failed_setup = 0;
    if ((res == false) || sim_state->simulation_should_kill_itself) {
        safe_log(temp_log, "ERROR! Something bad happened when reading the input file\n");
        failed_setup = 1;
    }

    if (sim_env->geometry == 0) {
        safe_log(temp_log, "ERROR! Structure type was not specified in input file\n");
        failed_setup = 1;
    }

    if (failed_setup) {
        safe_log(temp_log, "Saving log file for debugging\n");
        char temp_filename[256];
        snprintf(temp_filename, 256, "mesosim_temp_%d.log", (int)starttime);
        FILE *debug_file = fopen(temp_filename, "w+");
        if (debug_file == NULL) {
            fprintf(stderr, "Failed to save temporary log file: %s\n", strerror(errno));
            clean_and_error(errno);
        }
        write_backlog(temp_log, debug_file);
        fclose(debug_file);
        fclose(temp_log);
        clean_and_error(EXIT_FAILURE);
    }

    printf("Read file successfully\n");

    // put everything that was in temp_log into outFile
    write_backlog(temp_log, log_state->sim_log);
    fclose(temp_log);

    // finish_preprocessing();   //only called when deposition matters

    initialize_simulation(sim_state, sim_env, log_state);

    if (log_state->verbose)
        printf("Beginning simulation\n");

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

    // finalize everything

    time(&endtime);
    safe_log(log_state->sim_log, "Finished! Total time taken: %d seconds\n",
            (int)(endtime - starttime));

    clean_and_error(EXIT_SUCCESS);

    return 0;
}
