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

//use main function to run everything
int main(int argc, char* argv[]) {

    printf("main entered\n");
    // TODO: check if argc<=1 and print a usage blub and exit/return
    //start the time
    time(&starttime);

    //write to a temporary file until a logfile is identified
    // TODO: remove requirement of tmp folder
    char* temp_name = "tmp/tmp.log";
    temp_log = fopen(temp_name, "w");
    printf("Temp log: %s\n", temp_name);
    if (temp_log == NULL) {
        perror("Failed to make temp file");
        exit(errno);
    }
    fputs("MESOSIM 2024\n", temp_log);
    initialize_lattice_geometry(); //this gets overwritten by info from the input file

    simulation_type = -1; //TODO: need to define in globals!!
    // TODO: check length of argv with argc
    fprintf(temp_log, "Attempting to read in file %s\n", argv[1]);
    // get_input_file also initializes atom list
    if (get_input_file(argv[1]) == false) {
        fprintf(temp_log, "ERROR! Something bad happened when reading the input file\n");
        return 1;
    }

    if (simulation_should_kill_itself) {
        fprintf(temp_log, "ERROR! Something bad happened when reading the input file\n");
        return 1;
    }

    if (simulation_type == -1) {
        fprintf(temp_log, "ERROR! Structure type was not specified in input file\n");
        return 1;
    }

    printf("Read file successfully, finishing preprocessing now\n");
    //pre-process the file information and fill in the gaps with defaults

    //finish_preprocessing();   //only called when deposition matters

    general_simulation_initialization();

    do_initialize_simulation(simulation_type);


    if (strcmp(outFile, "") == 0) {
        //an out file name was not defined in input file, use starttime as filename
        sprintf(outFile, "%lld.out", starttime);
        fprintf(temp_log, "Log file name not defined, using \"%s\"", outFile);
    }

    sim_log_file = fopen(outFile, "w+");

    if (sim_log_file == NULL) {
        printf("ERROR! Could not open output file %s", outFile);
        return 1;
    }

    // put everything that was in temp_log into outFile
    write_backlog(temp_log, sim_log_file);
    fclose(temp_log);
    remove(temp_name);

    //print a lot of information to the log
    // TODO: move this to a function, esp since most of these are globals anyways
    fprintf(sim_log_file, "successfully read input file and preprocessed\n");
    fprintf(sim_log_file, "system size is %lf x %lf x %lf\n", ssx, ssy, ssz);

    switch (lattice_type) {
        case FCC:
            fprintf(sim_log_file, "crystal structure is FCC\n");
            break;
        case BCC:
            fprintf(sim_log_file, "crystal structure is BCC\n");
            break;
        case SC:
            fprintf(sim_log_file, "crystal structure is SC\n");
            break;
    }

    if (substrate_percent_a >= 1.) {
        //unary system
        fprintf(sim_log_file, "Initializing atom type %s\n", atom_names[0]);
        fprintf(sim_log_file, "Atom %s:\tsolubility %d;\tnearest neighbor energy [eV] %lf\n", atom_names[0], solubility[0], nnE[0]);
    }
    else if (substrate_percent_a + substrate_percent_b >= 1.-1e-4) {
        //binary system
        fprintf(sim_log_file, "Initializing atom types %s %s\n", atom_names[0], atom_names[1]);
        fprintf(sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf\n", atom_names[0], solubility[0], nnE[0], nnE[1]);
        fprintf(sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf\n", atom_names[1], solubility[1], nnE[1], nnE[3]);
    }
    else {
        //ternary system
        fprintf(sim_log_file, "Initializing atom types %s %s %s\n", atom_names[0], atom_names[1], atom_names[2]);
        fprintf(sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", atom_names[0], solubility[0], nnE[0], nnE[1], nnE[2]);
        fprintf(sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", atom_names[1], solubility[1], nnE[1], nnE[3], nnE[4]);
        fprintf(sim_log_file, "Atom %s:\tsolubility %d\tnearest neighbor energies [eV] %lf %lf %lf\n", atom_names[2], solubility[2], nnE[2], nnE[4], nnE[5]);
    }

    fprintf(sim_log_file, "Temperature is %lf K\n", temperature);

    if (overpotential_ramp_rate > 0.)
        fprintf(sim_log_file, "Potential sweep [eV/s] from %lf to %lf at %lf\n", overpotential, overpotential_ramp_rate, maxoverpotential);
    else
        fprintf(sim_log_file, "Potential constant [eV] at %lf\n", overpotential);

    if (analysis_type == REGULAR_TIME_INTERVALS)
        fprintf(sim_log_file, "Recording data at linear intervals [s] from %lf to %lf at %lf increments\n", time_interval_end, data_time_interval, run_time);
    else if (analysis_type == LOG_TIME_INTERVALS)
        fprintf(sim_log_file, "Recording data at log intervals [s] from %lf to %lf at %lf multiples\n", initial_logtime, logtime_multiplier, run_time);

    fprintf(sim_log_file, "Random seed is %ld\n", rand_seed);

    switch (simulation_type) {
        case SIMULATION_TYPE_FLAT_SHEET:
            fprintf(sim_log_file, "Initialized flat sheet with monolayer depth %d\n", sheet_thickness);
            break;
        case SIMULATION_TYPE_CLUSTER:
            fprintf(sim_log_file, "Initialized spherical cluster with radius %d\n", cluster_radius);
            break;
        case SIMULATION_TYPE_FROM_FILE:
            fprintf(sim_log_file, "Initialized user-defined structure with filename %s\n", atoms_filename);
            break;
    }

    fprintf(sim_log_file, "Atoms created, %d total\n", atom_cnt);
    fprintf(sim_log_file, "Beginning simulation\n");  
    
    //get the necessary file name prefix for xyz outputs
    strcpy(coordinate_log_prefix, outFile);
    char* lastdot = strrchr(coordinate_log_prefix, '.');
    lastdot[0] = '\0';


    //perform simulations
    unsigned long sim_error;
    sim_error = perform_simulation();

    if (sim_error != 0) {
        printf("ERROR! Something went wrong in the simulation\n");
        return 1;
    }

    //finalize everything

    //free the only malloced memory
    for (int i = transition_cnt; i > 0; --i)
        free(transition_arr[i-1]);
    for (int i = atom_cnt; i > 0; --i)
        free(atom_arr[i-1]);

    time(&endtime);
    fprintf(sim_log_file, "Finished! Total time taken: %d seconds\n", (int)(endtime-starttime));
    fclose(sim_log_file);
    return 0;
}
// initializes primitve_basis, ucell_params
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

	ssx = 1;
	ssy = 1;
	ssz = 1;

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