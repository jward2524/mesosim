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
    exit(0);
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
                fprintf(stderr, "Only one input file is supported: %s is the second unsupported argument\n", arg);
            }
            else {
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

    sim_state = calloc(1, sizeof(struct SimulationState));
    sim_env = calloc(1, sizeof(struct SimulationEnv));
    log_state = calloc(1, sizeof(struct LoggingState));

    set_state(sim_state, sim_env, log_state);

    // initialize_lattice_geometry(); //this gets overwritten by info from the input file
    sim_env->geometry = -1; // TODO: need to define in globals!!

    log_state->verbose = verbose_flag;

    // write to a temporary file until a logfile is identified
    //  not supported for msvcrt.dll [msys's mingw64]
    FILE *temp_log = tmpfile();
    if (temp_log == NULL) {
        perror("Failed to make temp file");
        clean_and_exit(errno);
    }
    if (log_state->verbose)
        printf("Temporary log created\n");
    fputs("MESOSIM 2024\n", temp_log);
    fprintf(temp_log, "Start time: %s\n", ctime(&starttime));
    fprintf(temp_log, "Attempting to read in file %s\n", input_filename);

    // TODO: clean up - move error handling into simulation_parameters_from_file; move all
    // initializers into one function in Sim_Aux simulation_parameters_from_file also initializes
    // atom list

    // pre-process the file information and fill in the gaps with defaults
    bool res = simulation_parameters_from_file(input_filename, sim_state, sim_env, log_state,
                                               temp_log, starttime);
    if (res == false) {
        fprintf(temp_log, "ERROR! Something bad happened when reading the input file\n");
        return 1;
    }

    if (sim_state->simulation_should_kill_itself) {
        fprintf(temp_log, "ERROR! Something bad happened when reading the input file\n");
        return 1;
    }

    if (sim_env->geometry == -1) {
        fprintf(temp_log, "ERROR! Structure type was not specified in input file\n");
        return 1;
    }

    printf("Read file successfully\n");

    // put everything that was in temp_log into outFile
    write_backlog(temp_log, log_state->sim_log_file);
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
        printf("ERROR! Something went wrong in the simulation\n");
        return 1;
    }

    // finalize everything

    time(&endtime);
    fprintf(log_state->sim_log_file, "Finished! Total time taken: %d seconds\n",
            (int)(endtime - starttime));

    clean_and_exit(0);

    return 0;
}
