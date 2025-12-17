#include "Mesosim.h"
#include "FileIO.h"
#include "Simulation.h"
#include "ErrorM.h"
#include "Initialization.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>

time_t starttime = 0;
time_t endtime = 0;

struct SimulationState *sim_state;
struct SimulationEnv *sim_env;
struct LoggingState *log_state;

void usage(){
    printf(
        "Usage: mesosim [FILE]\n"
        "Execute the mesosim KMC simulation program, with FILE as the input file\n"
        "\n"
        "Options:\n"
        "  -h, --help\tDisplay this message\n"
    );
    exit(0);
}

//use main function to run everything
int main(int argc, char* argv[]) {

    if ((argc <= 1) || (strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)){
        usage();
    }

    //start the time
    time(&starttime);
    
    sim_state = calloc(1, sizeof(struct SimulationState));
    sim_env = calloc(1, sizeof(struct SimulationEnv));
    log_state = calloc(1, sizeof(struct LoggingState));

    set_state(sim_state, sim_env, log_state);

    sim_env->zone_count_u = TTS;
    sim_env->zone_count_v = TTS;
    sim_env->zone_count_w = TTS;
    sim_env->dissolution = DISSOLUTION; // holy shit, when not =DISSOLUTION[=1], it breaks some shit heavy
    // initialize_lattice_geometry(); //this gets overwritten by info from the input file
    sim_env->simulation_type = -1; //TODO: need to define in globals!!
    
    //write to a temporary file until a logfile is identified
    // not supported for msvcrt.dll [msys's mingw64]
    FILE *temp_log = tmpfile();
    printf("Temporary log created\n");
    if (temp_log == NULL) {
        perror("Failed to make temp file");
        clean_and_exit(errno);
    }
    fputs("MESOSIM 2024\n", temp_log);
    fprintf(temp_log, "Start time: %lld\n", starttime);
    fprintf(temp_log, "Attempting to read in file %s\n", argv[1]);

    // TODO: clean up - move error handling into simulation_parameters_from_file; move all initializers into one function in Sim_Aux
    // simulation_parameters_from_file also initializes atom list
    if (simulation_parameters_from_file(argv[1], sim_state, sim_env, log_state, temp_log, starttime) == false) {
        fprintf(temp_log, "ERROR! Something bad happened when reading the input file\n");
        return 1;
    }

    if (sim_state->simulation_should_kill_itself) {
        fprintf(temp_log, "ERROR! Something bad happened when reading the input file\n");
        return 1;
    }

    if (sim_env->simulation_type == -1) {
        fprintf(temp_log, "ERROR! Structure type was not specified in input file\n");
        return 1;
    }

    printf("Read file successfully\n");
    //pre-process the file information and fill in the gaps with defaults

    // put everything that was in temp_log into outFile
    write_backlog(temp_log, log_state->sim_log_file);
    fclose(temp_log);

    //finish_preprocessing();   //only called when deposition matters
    
	// system geometry initialization

    // atom_cnt=0 for the initialization functions, so some of them end up doing nothing
    // bit shifts for periodic boundary conditions
	get_shifts(sim_env);
    
	set_primitive_basis(sim_env);
	set_default_orientation(sim_state->atom_arr, sim_state->atom_cnt, sim_env->lattice_type, sim_env->rmat, sim_env->primitive_basis); // supposedly was only for visualization
	get_system_normal(sim_env->primitive_basis);	// maybe only for visualization
    
    initialize_simulation_box(sim_env);
	
    // initialize zones - help figure out which atoms are next to which other atoms
	initialize_zones(sim_state->zone_arr, sim_env);
    
    // mallocs and sets to zero (or something else) simulation variables
    initialize_simulation_variables(sim_state, sim_env);
    
    // sets jump offsets for given crystal type
	initialize_neighbor_offsets(sim_env);
    initialize_initial_structure(sim_state, sim_env);
    check_system(sim_state, sim_env);

    // if (strcmp(outFile, "") == 0) {
    //     //an out file name was not defined in input file, use starttime as filename
    //     sprintf(outFile, "%lld.out", starttime);
    //     fprintf(temp_log, "Log file name not defined, using \"%s\"", outFile);
    // }

    // log_state->sim_log_file = fopen(outFile, "w+");

    // if (log_state->sim_log_file == NULL) {
    //     printf("ERROR! Could not open output file %s", outFile);
    //     return 1;
    // }

    // // put everything that was in temp_log into outFile
    // write_backlog(temp_log, log_state->sim_log_file);
    // fclose(temp_log);
    // remove(temp_name);

    input_logging(sim_state, sim_env, log_state);
    
    printf("Beginning simulation\n");
    //perform simulations
    unsigned long sim_error;
    sim_error = perform_simulation(sim_state, sim_env, log_state);

    if (sim_error != 0) {
        printf("ERROR! Something went wrong in the simulation\n");
        return 1;
    }

    //finalize everything

    time(&endtime);
    fprintf(log_state->sim_log_file, "Finished! Total time taken: %d seconds\n", (int)(endtime-starttime));
    
    clean_and_exit(0);

    return 0;
}
