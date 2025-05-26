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
    sim_env->dissolution = DISSOLUTION; // holy shit, when not =1, it breaks some shit heavy
    // Atom *arr[];
    // sim_state->atom_arr = &arr;
    // sim_env->max_neighbors = 12;
    // initialize_lattice_geometry(); //this gets overwritten by info from the input file

    sim_env->simulation_type = -1; //TODO: need to define in globals!!

    // TODO: multiple instances of program will be overwriting temp files
	//write to a temporary file until a logfile is identified
    char* temp_name = "temp.log";
    FILE* temp_log = fopen(temp_name, "w");
    printf("Temp log: %s\n", temp_name);
    if (temp_log == NULL) {
        perror("Failed to make temp file");
        exit(errno);
    }
    fputs("MESOSIM 2024\n", temp_log);
	fprintf(temp_log, "Start time: %lld", starttime);
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

    printf("Read file successfully, finishing preprocessing now\n");
    //pre-process the file information and fill in the gaps with defaults

    // put everything that was in temp_log into outFile
    write_backlog(temp_log, log_state->sim_log_file);
    fclose(temp_log);
    remove(temp_name);

    //finish_preprocessing();   //only called when deposition matters

    general_simulation_initialization(sim_state, sim_env, log_state);
    
    initialize_simulation_box(sim_env);
    // get_shifts()
    // initialize_zones()

    do_initialize_simulation(sim_state, sim_env);


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
    fprintf(log_state->sim_log_file, "system size is %lf x %lf x %lf\n", sim_env->ssx, sim_env->ssy, sim_env->ssz);

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

    if (sim_env->substrate_percent_a >= 1.) {
        //unary system
        fprintf(log_state->sim_log_file, "Initializing atom type %s\n", sim_env->atom_names[0]);
        fprintf(log_state->sim_log_file, "Atom %s:\tsim_env->solubility %d;\tnearest neighbor energy [eV] %lf\n", sim_env->atom_names[0], sim_env->solubility[0], sim_env->nnE[0]);
    }
    else if (sim_env->substrate_percent_a + sim_env->substrate_percent_b >= 1.-1e-4) {
        //binary system
        fprintf(log_state->sim_log_file, "Initializing atom types %s %s\n", sim_env->atom_names[0], sim_env->atom_names[1]);
        fprintf(log_state->sim_log_file, "Atom %s:\tsim_env->solubility %d\tnearest neighbor energies [eV] %lf %lf\n", sim_env->atom_names[0], sim_env->solubility[0], sim_env->nnE[0], sim_env->nnE[1]);
        fprintf(log_state->sim_log_file, "Atom %s:\tsim_env->solubility %d\tnearest neighbor energies [eV] %lf %lf\n", sim_env->atom_names[1], sim_env->solubility[1], sim_env->nnE[1], sim_env->nnE[3]);
    }
    else {
        //ternary system
        fprintf(log_state->sim_log_file, "Initializing atom types %s %s %s\n", sim_env->atom_names[0], sim_env->atom_names[1], sim_env->atom_names[2]);
        fprintf(log_state->sim_log_file, "Atom %s:\tsim_env->solubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", sim_env->atom_names[0], sim_env->solubility[0], sim_env->nnE[0], sim_env->nnE[1], sim_env->nnE[2]);
        fprintf(log_state->sim_log_file, "Atom %s:\tsim_env->solubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", sim_env->atom_names[1], sim_env->solubility[1], sim_env->nnE[1], sim_env->nnE[3], sim_env->nnE[4]);
        fprintf(log_state->sim_log_file, "Atom %s:\tsim_env->solubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", sim_env->atom_names[2], sim_env->solubility[2], sim_env->nnE[2], sim_env->nnE[4], sim_env->nnE[5]);
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
    fprintf(log_state->sim_log_file, "Beginning simulation\n");  
    
    //perform simulations
    unsigned long sim_error;
    sim_error = perform_simulation(sim_state, sim_env, log_state);

    if (sim_error != 0) {
        printf("ERROR! Something went wrong in the simulation\n");
        return 1;
    }

    //finalize everything

    //free the only malloced memory
    for (int i = sim_state->transition_cnt; i > 0; --i)
        free(sim_state->transition_arr[i-1]);
    for (int i = sim_state->atom_cnt; i > 0; --i)
        free(sim_state->atom_arr[i-1]);

    time(&endtime);
    fprintf(log_state->sim_log_file, "Finished! Total time taken: %d seconds\n", (int)(endtime-starttime));
    fclose(log_state->sim_log_file);
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

    double primitive_basis[3][3] = {
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.},
    };

	inver(primitive_basis, invert_primitive_basis);
	primitive_basis2ucell_params(primitive_basis, ucell_params);

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
