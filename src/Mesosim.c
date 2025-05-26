#include "stdafx.h"
#include "Defs.h"
#include "Geometry.h"
#include "Mesosim.h"
#include "Vector.h"
#include "Random.h"
#include <errno.h>
#include "FileIO.h"
#include "Simulation_Aux.h"
#include "Simulation.h"
#include "Atoms.h"
//notable absence of Windows.h!

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

    // TODO: multiple instances of program will be overwriting temp files
    //write to a temporary file until a logfile is identified
    char* temp_name = "temp.log";
    g_temp_log = fopen(temp_name, "w");
    printf("Temp log: %s\n", temp_name);
    if (g_temp_log == NULL) {
        perror("Failed to make temp file");
        exit(errno);
    }
    fputs("MESOSIM 2024\n", g_temp_log);
    initialize_lattice_geometry(); //this gets overwritten by info from the input file

    g_simulation_type = -1; //TODO: need to define in globals!!
    fprintf(g_temp_log, "Start time: %lld", starttime);
    fprintf(g_temp_log, "Attempting to read in file %s\n", argv[1]);
    // simulation_parameters_from_file also initializes atom list
    if (simulation_parameters_from_file(argv[1], sim_state) == false) {
        fprintf(g_temp_log, "ERROR! Something bad happened when reading the input file\n");
        return 1;
    }

    if (sim_state->simulation_should_kill_itself) {
        fprintf(g_temp_log, "ERROR! Something bad happened when reading the input file\n");
        return 1;
    }

    if (g_simulation_type == -1) {
        fprintf(g_temp_log, "ERROR! Structure type was not specified in input file\n");
        return 1;
    }

    printf("Read file successfully, finishing preprocessing now\n");
    //pre-process the file information and fill in the gaps with defaults

    //finish_preprocessing();   //only called when deposition matters

    general_simulation_initialization(sim_state);
    
    initialize_simulation_box(g_ssx, g_ssy, g_ssz);
    // get_shifts()
    // initialize_zones(, ss)

    do_initialize_simulation(g_simulation_type, sim_state);


    if (strcmp(g_outFile, "") == 0) {
        //an out file name was not defined in input file, use starttime as filename
        sprintf(g_outFile, "%lld.out", starttime);
        fprintf(g_temp_log, "Log file name not defined, using \"%s\"", g_outFile);
    }

    g_sim_log_file = fopen(g_outFile, "w+");

    if (g_sim_log_file == NULL) {
        printf("ERROR! Could not open output file %s", g_outFile);
        return 1;
    }

    // put everything that was in g_temp_log into g_outFile
    write_backlog(g_temp_log, g_sim_log_file);
    fclose(g_temp_log);
    remove(temp_name);

    //print a lot of information to the log
    // TODO: move this to a function, esp since most of these are globals anyways
    fprintf(g_sim_log_file, "successfully read input file and preprocessed\n");
    fprintf(g_sim_log_file, "system size is %lf x %lf x %lf\n", g_ssx, g_ssy, g_ssz);

    switch (g_lattice_type) {
        case FCC:
            fprintf(g_sim_log_file, "crystal structure is FCC\n");
            break;
        case BCC:
            fprintf(g_sim_log_file, "crystal structure is BCC\n");
            break;
        case SC:
            fprintf(g_sim_log_file, "crystal structure is SC\n");
            break;
    }

    if (g_substrate_percent_a >= 1.) {
        //unary system
        fprintf(g_sim_log_file, "Initializing atom type %s\n", g_atom_names[0]);
        fprintf(g_sim_log_file, "Atom %s:\tsolubility %d;\tnearest neighbor energy [eV] %lf\n", g_atom_names[0], g_solubility[0], g_nnE[0]);
    }
    else if (g_substrate_percent_a + g_substrate_percent_b >= 1.-1e-4) {
        //binary system
        fprintf(g_sim_log_file, "Initializing atom types %s %s\n", g_atom_names[0], g_atom_names[1]);
        fprintf(g_sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf\n", g_atom_names[0], g_solubility[0], g_nnE[0], g_nnE[1]);
        fprintf(g_sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf\n", g_atom_names[1], g_solubility[1], g_nnE[1], g_nnE[3]);
    }
    else {
        //ternary system
        fprintf(g_sim_log_file, "Initializing atom types %s %s %s\n", g_atom_names[0], g_atom_names[1], g_atom_names[2]);
        fprintf(g_sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", g_atom_names[0], g_solubility[0], g_nnE[0], g_nnE[1], g_nnE[2]);
        fprintf(g_sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", g_atom_names[1], g_solubility[1], g_nnE[1], g_nnE[3], g_nnE[4]);
        fprintf(g_sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", g_atom_names[2], g_solubility[2], g_nnE[2], g_nnE[4], g_nnE[5]);
    }

    fprintf(g_sim_log_file, "Temperature is %lf K\n", sim_state->temperature);

    if (g_overpotential_ramp_rate > 0.)
        fprintf(g_sim_log_file, "Potential sweep [eV/s] from %lf to %lf at %lf\n", sim_state->overpotential, g_overpotential_ramp_rate, g_max_overpotential);
    else
        fprintf(g_sim_log_file, "Potential constant [eV] at %lf\n", sim_state->overpotential);

    if (g_analysis_type == REGULAR_TIME_INTERVALS)
        fprintf(g_sim_log_file, "Recording data at linear intervals [s] from %lf to %lf at %lf increments\n", g_next_log_checkpoint, g_log_interval, sim_state->run_stime);
    else if (g_analysis_type == LN_TIME_INTERVALS)
        fprintf(g_sim_log_file, "Recording data at log intervals [s] from %lf to %lf at %lf multiples\n", g_next_log_checkpoint, g_log_interval, sim_state->run_stime);
    // TODO: fill out for other analysis_types

    fprintf(g_sim_log_file, "Random seed is %ld\n", rand_seed);

    switch (g_simulation_type) {
        case SIMULATION_TYPE_FLAT_SHEET:
            fprintf(g_sim_log_file, "Initialized flat sheet with monolayer depth %d\n", g_sheet_thickness);
            break;
        case SIMULATION_TYPE_CLUSTER:
            fprintf(g_sim_log_file, "Initialized spherical cluster with radius %d\n", g_cluster_radius);
            break;
        case SIMULATION_TYPE_FROM_FILE:
            fprintf(g_sim_log_file, "Initialized user-defined structure with filename %s\n", g_atoms_filename);
            break;
    }

    fprintf(g_sim_log_file, "Atoms created, %d total\n", g_atom_cnt);
    fprintf(g_sim_log_file, "Beginning simulation\n");  
    
    //get the necessary file name prefix for xyz outputs
    strcpy(g_coordinate_log_prefix, g_outFile);
    char* lastdot = strrchr(g_coordinate_log_prefix, '.');
    lastdot[0] = '\0';


    //perform simulations
    unsigned long sim_error;
    sim_error = perform_simulation(sim_state);

    if (sim_error != 0) {
        printf("ERROR! Something went wrong in the simulation\n");
        return 1;
    }

    //finalize everything

    //free the only malloced memory
    for (int i = g_transition_cnt; i > 0; --i)
        free(g_transition_arr[i-1]);
    for (int i = g_atom_cnt; i > 0; --i)
        free(g_atom_arr[i-1]);

    time(&endtime);
    fprintf(g_sim_log_file, "Finished! Total time taken: %d seconds\n", (int)(endtime-starttime));
    fclose(g_sim_log_file);
    return 0;
}
// initializes primitve_basis, ucell_params, ss*
void initialize_lattice_geometry(void)
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

	g_ssx = 1;
	g_ssy = 1;
	g_ssz = 1;

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
