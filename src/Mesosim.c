#include "Mesosim.h"
#include "FileIO.h"
#include "Random.h"
#include "Simulation_Aux.h"
#include "Simulation.h"
#include "Atoms.h"
#include "Vector.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>

time_t starttime = 0;
time_t endtime = 0;

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
// TODO: create full/safe exit function that closes files somewhere (closes files + deletes temp)

//use main function to run everything
int main(int argc, char* argv[]) {

    if ((argc <= 1) || (strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)){
        usage();
    }

    //start the time
    time(&starttime);
    
    struct SimulationState *sim_state = calloc(1, sizeof(struct SimulationState));
    struct SimulationEnv *sim_env = calloc(1, sizeof(struct SimulationEnv));
    struct LoggingState *log_state = calloc(1, sizeof(struct LoggingState));

    sim_env->zix = TTS;
    sim_env->ziy = TTS;
    sim_env->ziz = TTS;
    sim_env->dissolution = DISSOLUTION; // holy shit, when not =DISSOLUTION[=1], it breaks some shit heavy
    // initialize_lattice_geometry(); //this gets overwritten by info from the input file
    sim_env->simulation_type = -1; //TODO: need to define in globals!!
    
    //write to a temporary file until a logfile is identified
    // not supported for msvcrt.dll [msys's mingw64]
    FILE *temp_log = tmpfile();
    printf("Temporary log created\n");
    if (temp_log == NULL) {
        perror("Failed to make temp file");
        exit(errno);
    }
    fputs("MESOSIM 2024\n", temp_log);
    fprintf(temp_log, "Start time: %lld\n", starttime);
    fprintf(temp_log, "Attempting to read in file %s\n", argv[1]);

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

    general_simulation_initialization(sim_state, sim_env, log_state);
    
    initialize_simulation_box(sim_env);
    // get_shifts()
    // initialize_zones()

    initialize_initial_structure(sim_state, sim_env);


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

    //print a lot of information to the log
    // TODO: move this to a function, esp since most of these are globals anyways
    fprintf(log_state->sim_log_file, "successfully read input file and preprocessed\n");
    fprintf(log_state->sim_log_file, "system size is %d x %d x %d\n", sim_env->ssx, sim_env->ssy, sim_env->ssz);

    switch (sim_env->lattice_type) {
        case FCC:
            fprintf(log_state->sim_log_file, "crystal structure is FCC\n");
            break;
        case BCC:
            fprintf(log_state->sim_log_file, "crystal structure is BCC\n");
            break;
        case SC:
            fprintf(log_state->sim_log_file, "crystal structure is SC\n");
            break;
    }

    fprintf(log_state->sim_log_file, "Initializing atom types: ");
    for (int i = 0; i < sim_env->num_elements; i++)
    {
        fprintf(log_state->sim_log_file, "%s ", sim_env->atom_names[i]);
    }
    fprintf(log_state->sim_log_file, "\nComposition: ");
    for (int i = 0; i < sim_env->num_elements; i++)
    {
        fprintf(log_state->sim_log_file, "%lf ", sim_env->substrate_composition[i]);
    }

    fprintf(log_state->sim_log_file, "\nSolutility: ");
    for (int i = 0; i < sim_env->num_elements; i++)
    {
        fprintf(log_state->sim_log_file, "%s ", sim_env->is_soluble[i] ? "false" : "true");
    }

    fprintf(log_state->sim_log_file, "\nBond energies\n");
    int bond_idx, env_idx;
    for (int i = 0; i < sim_env->num_elements; i++)
    {
        for (int j = i; j < sim_env->num_elements; j++)
        {
            bond_idx = get_bond_index(i, j, sim_env);
            env_idx = get_env_index(1, bond_idx, sim_env);
            fprintf(log_state->sim_log_file, "%s-%s: %lf\n", sim_env->atom_names[i], sim_env->atom_names[j], sim_env->nn_energy[env_idx]);
        }
    }

    fprintf(log_state->sim_log_file, "Temperature is %lf K\n", sim_env->temperature);

    if (sim_env->overpotential_ramp_rate > 0.)
        fprintf(log_state->sim_log_file, "Potential sweep [eV/s] from %lf to %lf at %lf\n", sim_env->overpotential, sim_env->overpotential_ramp_rate, sim_env->max_overpotential);
    else
        fprintf(log_state->sim_log_file, "Potential constant [eV] at %lf\n", sim_env->overpotential);

    if (log_state->analysis_type == REGULAR_TIME_INTERVALS)
        fprintf(log_state->sim_log_file, "Recording data at linear intervals [s] from %lf to %lf at %lf increments\n", log_state->next_log_checkpoint, log_state->log_interval, sim_state->run_stime);
    else if (log_state->analysis_type == LN_TIME_INTERVALS)
        fprintf(log_state->sim_log_file, "Recording data at log intervals [s] from %lf to %lf at %lf multiples\n", log_state->next_log_checkpoint, log_state->log_interval, sim_state->run_stime);
    // TODO: fill out for other analysis_types

    fprintf(log_state->sim_log_file, "Random seed is %ld\n", rand_seed);

    switch (sim_env->simulation_type) {
        case SIMULATION_TYPE_FLAT_SHEET:
            fprintf(log_state->sim_log_file, "Initialized flat sheet with monolayer depth %d\n", sim_env->sheet_thickness);
            break;
        case SIMULATION_TYPE_CLUSTER:
            fprintf(log_state->sim_log_file, "Initialized spherical cluster with radius %d\n", sim_env->cluster_radius);
            break;
        case SIMULATION_TYPE_FROM_FILE:
            fprintf(log_state->sim_log_file, "Initialized user-defined structure with filename %s\n", sim_env->atoms_filename);
            break;
    }

    fprintf(log_state->sim_log_file, "Atoms created, %d total\n", sim_state->atom_cnt);
    
    printf("Beginning simulation\n");
    //perform simulations
    unsigned long sim_error;
    sim_error = perform_simulation(sim_state, sim_env, log_state);

    if (sim_error != 0) {
        printf("ERROR! Something went wrong in the simulation\n");
        return 1;
    }

    //finalize everything

    //free the malloc'ed memory
    for (int i = sim_state->transition_cnt; i > 0; --i)
    {
        free(sim_state->transition_arr[i-1]);
        sim_state->transition_arr[i-1] = NULL;
    }
    for (int i = sim_state->atom_cnt; i > 0; --i)
    {
        free(sim_state->atom_arr[i-1]);
        sim_state->atom_arr[i-1] = NULL;
    }

    
    free(sim_state->atom_arr);
    free(sim_state->rate_arr);
    free(sim_state->transition_arr);
    sim_state->atom_arr = NULL;
    sim_state->rate_arr = NULL;
    sim_state->transition_arr = NULL;

    free(sim_state);
    free(sim_env);
    sim_state = NULL;
    sim_env = NULL;
    
    time(&endtime);
    fprintf(log_state->sim_log_file, "Finished! Total time taken: %d seconds\n", (int)(endtime-starttime));
    fclose(log_state->sim_log_file);

    free(log_state);
    log_state = NULL;

    return 0;
}
// initializes primitve_basis, ucell_params, ss*
void initialize_lattice_geometry(struct SimulationEnv* sim_env)
{
	// Initializes the generic lattice geometry to be simple cubic (i.e., a=1, b=1, c=1, alpha = 90, beta = 90, gamma = 90)

	// int i,j;

	// for (i=0;i<3;++i)
	// 	for (j=0;j<3;++j)
	// 		if (i == j) primitive_basis[i][j] = 1.0; else primitive_basis[i][j] = 0.0;

    double pb[3][3] = {
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.},
    };
    memcpy(sim_env->primitive_basis, pb, 3*3*sizeof(double));

	inver(sim_env->primitive_basis, sim_env->invert_primitive_basis);
	primitive_basis2ucell_params(sim_env->primitive_basis, sim_env->ucell_params);

	sim_env->ssx = 1;
	sim_env->ssy = 1;
	sim_env->ssz = 1;

	return;
}

void write_backlog(FILE* tempFile, FILE* logFile)
{
    //transfers everything from the temporary log file to a permanent log
    rewind(tempFile);
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), tempFile) != NULL) {
        fputs(buffer, logFile);
    }
}
