#include "FileIO.h"
#include "Random.h"
#include "Atoms.h"
#include "Mesosim.h"
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>

const int BUFFER_SIZE = 200;
const int ARR_BUFFER_SIZE = 20;
char outFile[260] = ""; //MAX_PATH variable Windows related, default 260

int fact(int n);
void calloc_nnE(struct SimulationEnv* se);

bool fopen_error(char* filename, FILE* file, char* base_msg){
	if (file == NULL) 
	{
		fprintf(stderr, "%s%s: %s\n", base_msg, filename, strerror(errno));
        clean_and_exit(errno);
    }
	return true;
}

bool simulation_parameters_from_file(char *filename, struct SimulationState *ss, struct SimulationEnv *se, struct LoggingState *ls, FILE *temp_log, time_t starttime)
{	
	FILE* input_file = fopen(filename, "r");
	fopen_error(filename, input_file, "Failed to open input file, ");
	
	char extension[4] = "";
	char* file_ender = strrchr(filename, '.'); // everything after the final '.' in the filename
	// [ ]: what is happening here?
	// TODO: when restarting a simulation, input file will be a previous log file?

	bool ret = false;
	if (file_ender == NULL) // when is this null?
	{
		// if (sim_log_file == NULL) // when is this not null?
		fprintf(temp_log, "ERROR! extension not found in file: %s\n", filename); //for reading arguments
		// else
		// fprintf(sim_log_file, "ERROR! extension not found in file: %s\n", filename); //should only happen when restarting/checkpointing
	}
	else {
		file_ender++;
		strncpy(extension, file_ender, 3); // copy only extension into extension array
		fprintf(temp_log, "Reading input from .%s file, %s\n", extension, filename);
	}

	if (strncmp(extension, "xyz", 3) == 0)
	{
		// open simple x,y,z,type coordinate file
		ret = process_xyz_file(temp_log, input_file, ss, se, ls);
	}
	else if (strncmp(extension, "kmc", 3) == 0)
	{
		// open kmc type file
		ret = process_kmc_file(temp_log, input_file, ss, se, ls);
	}
	else if (strncmp(extension, "in", 2) == 0)
	{
		//open and process parameter input file
		ret = process_in_file(temp_log, input_file, ss, se, ls);
    }
	else if (strncmp(extension, "kmx", 3) == 0)
	{
		//open and process new kmc input type that removes fluff (does this need to happen)
		ret = process_kmx_file(temp_log, input_file, ss, se, ls);
    }
	else
	{
		// if (sim_log_file == NULL)
		fprintf(temp_log, "ERROR! file extension not recognized: .%s\n", extension); //for reading arguments
		// else
			// fprintf(sim_log_file, "ERROR! file extension not recognized: .%s\n", extension); //should only happen when restarting/checkpointing
		ret = false;
    }
	
	if (strcmp(outFile, "") == 0) {
		//an out file name was not defined in input file, use starttime as filename
        sprintf(outFile, "%lld.out", starttime);
        fprintf(temp_log, "Log file name not defined, using \"%s\"", outFile);
    }
	
    ls->sim_log_file = fopen(outFile, "w+");
    ret = ret && fopen_error(filename, ls->sim_log_file, "Failed to open log file, ");

	//get the file name prefix for xyz outputs (outFile without the [.out] extension)
	strcpy(ls->position_log_prefix, outFile);
	char* lastdot = strrchr(ls->position_log_prefix, '.');
	lastdot[0] = '\0';
	return ret;
}


/*******************************************************************************
*******************************************************************************/

bool process_in_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls) {
	char parameter_line[BUFFER_SIZE];
	int errnum;
	// ENHANCE: line length should be a const that is used to pull lines and create buffer sizes
	while (fgets(parameter_line, BUFFER_SIZE, input_file) != NULL) {
		if (strncmp(parameter_line, "restart", 7) == 0) {
			//TODO: restart the simulation from a log file and don't do the rest of the loop
			// restart also needs an atom position file
		}
		errnum = parse_input(parameter_line, temp_log, ss, se, ls);
		if (errnum != SUCCESS)
		{
			//should write to the temp
			fprintf(temp_log, "ERROR! Had issue reading the following line: \"%s\"\n", parameter_line);
			fclose(input_file);
			return false;
		}
	}
	fclose(input_file);
	se->max_rates = (unsigned long long int) pow(se->max_neighbors + se->dissolution, se->num_nn_types);
	return true;
}

/*******************************************************************************
*******************************************************************************/

int parse_input(char* line, FILE* temp_log, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls)
{
	//TODO: use strtok?
	//printf("Trying to parse this line! \"%s\"\n", line);
	char* ptr = line; // line from file
	char cmd[BUFFER_SIZE]; // command - first word in line
	char params[BUFFER_SIZE]; // parameters parsed from line

	// line start with command word and is followed by parameters
	// split command word from parameters
	while (ptr[0] != '\n' && ptr[0] != '\0' && ptr[0] != ' ' && ptr[0] != '\t')
		++ptr; //find the first whitespace on the line (or the end of the line)

	strcpy(params, ptr);

	if (params[0] == '\0' || params[0] == '\n')
		return NOT_ENOUGH_PARAMS; //there are no parameters, at end of line


	// remove whitespace characters from left side of params
	int i = 0;
	while (params[i] != '\0' && (params[i] == ' ' || params[i] == '\t'))
		++i; //get rid of whitespace on left
	
	int k = 0;
	for (int j = i; params[j] != '\0'; ++j)
	{
		params[k] = params[j];
		++k;
	}
	params[k] = '\0'; // null-terminate the effective string

	if (params[0] == '\0' || params[0] == '\n')
		return NOT_ENOUGH_PARAMS; //there are no parameters, end of line
	
	sscanf(line, "%s", cmd); // puts `original` line into cmd

	//now handle individual keywords - check which one is at beginning of line=cmd
	int argsread; //check to see if everything got read correctly
	if (strncmp(cmd, "systemsize", 10) == 0) {
		// set the system size using params
		if ((argsread = sscanf(params, "%d %d %d", &se->system_size_x, &se->system_size_y, &se->system_size_z)) != 3)
		{
			fprintf(temp_log, "ERROR! Could not correctly read system size parameters %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
		se->max_atoms = se->system_size_x * se->system_size_y * se->system_size_z;
	}
	else if (strncmp(cmd, "temp", 4) == 0) {
		// set the system temperature
		if ((argsread = sscanf(params, "%lf", &ss->temperature)) != 1)
		{
			fprintf(temp_log, "ERROR! Could not correctly read temperature parameter %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
    }
	else if (strncmp(cmd, "seed", 4) == 0) {
		// set the random seed
		if (strncmp(params, "random", 6) == 0) {
			//rand_seed should be based on time
			time_t seedtime;
			time(&seedtime);
			rand_seed = (long int)seedtime;
			fprintf(temp_log, "Using random time seed %ld\n", rand_seed);
		}
		else if (strncmp(params, "default", 7) == 0) {
			rand_seed = DEFAULT_SEED;
			fprintf(temp_log, "Using default time seed %ld\n", rand_seed);
		}
		else if ((argsread = sscanf(params, "%ld", &rand_seed)) != 1) //read in a long int
		{
			fprintf(temp_log, "ERROR! Could not correctly read random seed parameter %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "potential", 9) == 0) {
		//set the potential sweep
		argsread = sscanf(params, "%lf %lf %lf", &se->initial_overpotential, &se->overpotential_ramp_rate, &se->max_overpotential);
		if (argsread == 1) {
			//constant overpotential
			se->overpotential_ramp_rate = 0.0;
			se->max_overpotential = se->initial_overpotential;
			fprintf(temp_log, "Using constant potential %lf\n", se->initial_overpotential);
		}
		else if (argsread != 3) //
		{
			fprintf(temp_log, "ERROR! Could not correctly read potential sweep parameters %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd+1, "nne", 3) == 0) {
		int nn_level;

		if ((se->num_nn_levels == 0) || (se->num_bond_types == 0))
		{
			fprintf(
				temp_log, 
				"ERROR! Number of nearest neighbor levels and number of elements \
				needs to be set before defining nearest neighbor energies\n"
			);
			return FILE_COMMAND_IGNORED;
		}

		int ret = sscanf(cmd, "%dnne", &nn_level);

		if (ret == 0)
		{
			fprintf(temp_log, "ERROR! Expected nne command format of [level]nne, recieved %s\n", cmd);
			return FILE_COMMAND_IGNORED;
		}
		
		if (nn_level > se->num_nn_levels)
		{
			fprintf(
				temp_log, 
				"ERROR! NN energy provided for higher level than stated: \
				%d stated nn levels, energies for %d level provided\n", 
				se->num_nn_levels, nn_level
			);
			return FILE_COMMAND_IGNORED;
		}

		// expect num_bond_types numbers
		if (se->nn_energy == NULL)
			calloc_nnE(se);
		int bond_index = get_env_index(nn_level, 0, se);
		int count = 0;

		int len = strlen(params); // BUFFER_SIZE?
		char tok_params[len];
		snprintf(tok_params, len, "%s", params);

		char* token = strtok(tok_params, " \t");
		while (token)
		{
			sscanf(token, "%lf", (se->nn_energy)+bond_index);
			token = strtok(NULL, " \t");
			bond_index++;
			count++;
		}

		if (count != se->num_bond_types)
		{
			fprintf(temp_log, "ERROR! Expected %d bond energies, recieved %d\n", se->num_bond_types, count);
			return FILE_COMMAND_IGNORED;
		}

	}
	else if (strncmp(cmd, "nnlevels", 8) == 0) 
	{
		argsread = sscanf(params, "%d", &(se->num_nn_levels));
		// if (se->num_bond_types != 0)
		// 	calloc_nnE(se);
	}
	else if (strncmp(cmd, "datalog", 4) == 0) { // ENHANCE: linear list and ln list do the same thing - improve semantics to eliminate this duplicity
		// set time increments for data logging
		int cursor;
		if (strncmp(params, "linear", 6) == 0) // linear data recording
		{ // ENHANCE: use cursor variable to step through params
			cursor = strlen("linear") + 1;
			int type = parse_datalog_params(params, cursor, ls, temp_log);
			ls->analysis_type = (type > 0) ? REGULAR_TIME_INTERVALS : TIME_LIST;
		}
		else if (strncmp(params, "ln", 2) == 0) // logarithmic data recording
		{
			cursor = strlen("ln") + 1;
			int type = parse_datalog_params(params, cursor, ls, temp_log);
			ls->analysis_type = (type > 0) ? LN_TIME_INTERVALS : TIME_LIST;
		}
		else if (strncmp(params, "iteration", 9) == 0) {
			cursor = strlen("iteration") + 1;
			int type = parse_datalog_params(params, cursor, ls, temp_log);
			ls->analysis_type = (type > 0) ? ITERATION_INTERVALS : ITERATION_LIST;
		}
		else {
			fprintf(temp_log, "ERROR! Unknown argument in 'datalog' command: %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	/*else if (strncmp(cmd, "runs", 4) == 0) {
		//set number of runs
		if ((argsread = sscanf(params, "%d", &number_of_simulation_runs)) != 1) {
			printf("ERROR! Could not correctly read number of simulation runs %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}*/
	else if (strncmp(cmd, "struct", 6) == 0) {
		//set crystal structure
		char structtype[3];
		if ((argsread = sscanf(params, "%s", structtype)) != 1) {
			fprintf(temp_log, "ERROR! Could not correctly read structure type %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
		if (strncmp(structtype, "FCC", 3) == 0 || strncmp(structtype, "1", 1) == 0)
			se->lattice_type = FCC;
		else if (strncmp(structtype, "BCC", 3) == 0 || strncmp(structtype, "2", 1) == 0)
			se->lattice_type = BCC;
		else if (strncmp(structtype, "SC", 2) == 0 || strncmp(structtype, "3", 1) == 0)
			se->lattice_type = SC;
		/*else if (strncmp(structtype, "DIA", 3) == 0 || strncmp(structtype, "4", 1) == 0)
			lattice_type = DIAMOND;*/ 
		else {
			fprintf(temp_log, "ERROR! Structure type %s not valid\n", structtype);
			return FILE_COMMAND_IGNORED;
		}
		se->max_neighbors = MAXIMUM_NUMBER_OF_NEIGHBORS;
	}
	else if (strncmp(cmd, "output", 6) == 0) {
		//set file name for log output
		if ((argsread = sscanf(params, "%s", outFile)) != 1) {
			fprintf(temp_log, "ERROR! Could not correctly read output file name %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "geometry", 8) == 0) {
		//initialize the atoms! the options are either flat sheet, spherical cluster, or file input
		if (strncmp(params, "sheet", 5) == 0) {
			se->simulation_type = SIMULATION_TYPE_FLAT_SHEET;
			if ((argsread = sscanf(params, "%*s %d", &se->sheet_thickness)) != 1) {
				fprintf(temp_log, "ERROR! Could not correctly read sheet thickness\n");
				return FILE_COMMAND_IGNORED;
			}
		}
		else if (strncmp(params, "cluster", 6) == 0) {
			se->simulation_type = SIMULATION_TYPE_CLUSTER;
			if ((argsread = sscanf(params, "%*s %d", &se->cluster_radius)) != 1) {
				fprintf(temp_log, "ERROR! Could not correctly read cluster radius\n");
				return FILE_COMMAND_IGNORED;
			}
		}
		else if (strncmp(params, "file", 4) == 0) {
			se->simulation_type = SIMULATION_TYPE_FROM_FILE;
			if ((argsread = sscanf(params, "%*s %s", se->atoms_filename)) != 1) {
				fprintf(temp_log, "ERROR! Could not correctly read file name for atoms\n");
				return FILE_COMMAND_IGNORED;
			}
		}
		else {
			fprintf(temp_log, "ERROR! Could not recognize geometry type %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "atomtype", 8) == 0) {
		//this determines which elements are which types of atoms

		char *types[ARR_BUFFER_SIZE];
		int len = strlen(params); // BUFFER_SIZE?
		char tok_params[len];
		snprintf(tok_params, len, "%s", params);
		char* token = strtok(tok_params, " \t");
		int count = 0;
		while (token)
		{
			types[count] = (char *)malloc(BUFFER_SIZE * sizeof(char)); // TODO: free // ENHANCE: 3 -> BUFFER_SIZE??
			sscanf(token, "%s", types[count]);
			token = strtok(NULL, " \t");
			count++;
		}

		if (count == 0){
			fprintf(temp_log, "ERROR! Couldn't read any atom type names %s\n", line);
			return FILE_COMMAND_IGNORED;
		}

		se->num_elements = count;
		se->num_bond_types = get_num_bond_types(se->num_elements);
		
		se->atom_names = (char **)calloc(count, sizeof(char*)); // TODO: free
		if (se->atom_names == NULL)
		{
			fprintf(stderr, "Couldn't allocate memory for atom names: %s", strerror(errno));
			fprintf(temp_log, "Couldn't allocate memory for atom names: %s", strerror(errno));
        	clean_and_exit(errno);
		}
		memcpy(se->atom_names, types, count * sizeof(char*));
		// if (se->num_nn_levels != 0)
		// 	calloc_nnE(se);

	}
	else if (strncmp(cmd, "dissolution", 11) == 0) {
		//this determines which atoms dissolve

		bool is_soluble[ARR_BUFFER_SIZE];
		int len = strlen(params); // BUFFER_SIZE?
		char tok_params[len];
		snprintf(tok_params, len, "%s", params);
		char* token = strtok(tok_params, " \t");
		int count = 0;
		while (token)
		{
			char buf[BUFFER_SIZE];
			sscanf(token, "%s", buf);
			int b = parse_boolean(buf);
			if (b < 0)
			{
				fprintf(temp_log, "ERROR! Could not correctly read solubility %s\n", buf);
				return FILE_COMMAND_IGNORED;
			}
			is_soluble[count] = (bool) b;
			token = strtok(NULL, " \t");
			count++;
		}

		if (count == 0)
		{
			fprintf(temp_log, "ERROR! Could not read any solubilities %s\n", line);
			return FILE_COMMAND_IGNORED;
		}

		if ((se->num_elements == 0) || (se->num_elements != count))
		{
			fprintf(temp_log, "ERROR! More values provided (%d) than number of elements %d\n", count, se->num_elements);
			return FILE_COMMAND_IGNORED;
		}

		int size = count * sizeof(bool);
		se->is_soluble = (bool *)malloc(size);
		if (se->is_soluble == NULL)
		{
			fprintf(stderr, "Couldn't allocate memory for solubilities: %s", strerror(errno));
			fprintf(temp_log, "Couldn't allocate memory for solubilities: %s", strerror(errno));
        	clean_and_exit(errno);
		}
		memcpy(se->is_soluble, is_soluble, size);
	}
	else if (strncmp(cmd, "composition", 11) == 0) {

		double comp[ARR_BUFFER_SIZE];
		int len = strlen(params); // BUFFER_SIZE?
		char tok_params[len];
		snprintf(tok_params, len, "%s", params);
		char* token = strtok(tok_params, " \t");
		int count = 0;
		double tot = 0;
		double c;
		while (token)
		{
			sscanf(token, "%lf", &c);
			comp[count] = c;
			tot += c;

			token = strtok(NULL, " \t");
			count++;
		}

		if (count == 0)
		{
			fprintf(temp_log, "ERROR! Could not read any compositions %s\n", line);
			return FILE_COMMAND_IGNORED;
		}

		if (count != se->num_elements)
		{
			fprintf(temp_log, "ERROR! More values provided (%d) than number of elements %d\n", count, se->num_elements);
			return FILE_COMMAND_IGNORED;
		}

		if (fabs(tot - 1) > 1e-10)
		{
			fprintf(temp_log, "ERROR! Compositions must add up to 1 - current: %lf\n", tot);
			return FILE_COMMAND_IGNORED;
		}

		int size = count * sizeof(double);
		se->substrate_composition = (double *)malloc(size);
		if (se->substrate_composition == NULL)
		{
			fprintf(stderr, "Couldn't allocate memory for compositions: %s", strerror(errno));
			fprintf(temp_log, "Couldn't allocate memory for compositions: %s", strerror(errno));
        	clean_and_exit(errno);
		}
		for (int i = 0; i < count; i++)
			se->substrate_composition[i] = comp[i];
	}
	else if (strncmp(cmd, "run", 3) == 0) {
		int cursor;
		if (strncmp(params, "time", 4) == 0) {
			cursor = 5;
			ss->sim_end_type = SIM_END_BY_STIME;
		}
		else if (strncmp(params, "iteration", 9) == 0) {
			cursor = 10;
			ss->sim_end_type = SIM_END_BY_ITERATIONS;
		}
		else {
			fprintf(temp_log, "ERROR! Unknown argument in 'run' command: %s\n", params);
			return FILE_COMMAND_IGNORED;
		}

		if (ss->sim_end_type == SIM_END_BY_STIME) {
			ss->run_stime = strtod(&params[cursor], NULL);
		}
		else if (ss->sim_end_type == SIM_END_BY_ITERATIONS) {
			ss->final_iteration = strtol(&params[cursor], NULL, 10);
		}

	}
	else {
		fprintf(temp_log, "ERROR! keyword %s not recognized\n", cmd);
		return FILE_COMMAND_IGNORED;
	}
	return SUCCESS;
}

int parse_boolean(char *str) {
	if (strncmp(str, "true", 4) == 0 || strncmp(str, "True", 4) == 0 || strncmp(str, "TRUE", 4) == 0 || strncmp(str, "T", 1) == 0 || strncmp(str, "1", 1) == 0)
		return 1;
	else if (strncmp(str, "false", 5) == 0 || strncmp(str, "False", 5) == 0 || strncmp(str, "FALSE", 5) == 0 || strncmp(str, "F", 1) == 0 || strncmp(str, "0", 1) == 0)
		return 0;
	else
		return -1;
}

void parse_log_list(char* input_str, double* list, int* len){
	char delim[] = " ";
	char* token = strtok(input_str, delim);
	while (token) {
		list[*len] = strtod(token, NULL);
		(*len)++;
		token = strtok(NULL, delim);
	}
}

// within datalog command, parsing interval and list keywords
// counter is 1 for simulation time, -1 for iterations (as per the macros)
// returns 1 if intervals, -1 if list
// ENHANCE: return values aren't conventional - make more conventional
int parse_datalog_params(char* params, int cursor, struct LoggingState* ls, FILE* temp_log){
	int argsread;
	if (strncmp(&params[cursor], "interval", 8) == 0) {
		argsread = sscanf(params + cursor + 9, "%lf %lf", &ls->next_log_checkpoint, &ls->log_interval);
		if (argsread < 2){
			ls->log_interval = ls->next_log_checkpoint;
		}
		return 1;
	}
	else if (strncmp(params + cursor, "list", 4) == 0){
		ls->log_list = (double*) malloc(ARR_BUFFER_SIZE * sizeof(double));
		if (ls->log_list == NULL)
		{
			fprintf(stderr, "Couldn't allocate memory for log list: %s", strerror(errno));
			fprintf(temp_log, "Couldn't allocate memory for log list: %s", strerror(errno));
        	clean_and_exit(errno);
		}
		parse_log_list(params + cursor + 5, ls->log_list, &ls->log_list_len);
		return -1;
	}
	else {
		fprintf(temp_log, "ERROR! Unknown argument in 'datalog' command: %s\n", params);
		return FILE_COMMAND_IGNORED;
	}
}

// allocates space for nnE array
void calloc_nnE(struct SimulationEnv* se)
{
	se->num_nn_types = se->num_nn_levels * se->num_bond_types;
	se->nn_energy = (double *)calloc(se->num_nn_types, sizeof(double));
}

// gets the index in atom_env, nnE arrays (nearest_neighbor - bond_type combo)
int get_env_index(int nn, int bond_idx, struct SimulationEnv* se)
{
	return (nn-1) * se->num_bond_types + bond_idx;
}

// int get_env_index_types(int nn, int a, int b, struct SimulationEnv* se)
// {
// 	return get_env_index(nn, get_bond_index(a, b, se), se);
// }

// calculate the number of bond types
int get_num_bond_types(int num_elements)
{
	return fact(num_elements + 2 - 1) / (fact(2) * fact(num_elements-1));
}

// returns the factorial of n
int fact(int n)
{
	switch(n)
	{
		case 0:
			return 1;
		case 1:
			return 1;
		case 2:
			return 2;
		case 3:
			return 6;
		case 4:
			return 24;
		case 5:
			return 120;
		case 6:
			return 720;
		case 7:
			return 5040;
		default:
			fprintf(stderr, "Factorial is too large (max n is 7): %d", n);
			return -1;
	}
}
/*******************************************************************************
*******************************************************************************/

bool process_kmc_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls)
{
	int newnat;

	int i,j,k;
	int x,y,z;

	// system and zone size
	se->system_size_x = DSIMSIZE; //is this always true?????
	se->system_size_y = DSIMSIZE;
	se->system_size_z = DSIMSIZE;
	
	// zones in x, y, z
	se->zone_count_u = TTS;
	se->zone_count_v = TTS;
	se->zone_count_w = TTS;	

	//general_simulation_initialization(); //happens later

	//read in the lattice and rotation matrices
	fscanf(input_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&se->primitive_basis[0][0], &se->primitive_basis[0][1], &se->primitive_basis[0][2], 
		&se->primitive_basis[1][0], &se->primitive_basis[1][1], &se->primitive_basis[1][2], 
		&se->primitive_basis[2][0], &se->primitive_basis[2][1], &se->primitive_basis[2][2]); 
		
	fscanf(input_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&se->rmat[0][0], &se->rmat[0][1], &se->rmat[0][2], 
		&se->rmat[1][0], &se->rmat[1][1], &se->rmat[1][2], 
		&se->rmat[2][0], &se->rmat[2][1], &se->rmat[2][2]); 

	fscanf(input_file, "%d %d %d", &se->system_size_x, &se->system_size_y, &se->system_size_z);

	fscanf(input_file, "%d", &newnat);

	fprintf(temp_log, "system size %d %d %d, number of atoms %d\n", se->system_size_x, se->system_size_y, se->system_size_z, newnat);
	
	int tempint;
	double tempdouble[3][3];
	Atom temp_atom;

	for (i=0;i<newnat;++i)
	{
		fscanf(input_file, "%s\t", temp_atom.name);

		fscanf(input_file, 
			"%c\t\
			%lf\t%lf\t%lf\t\
			%d\t%d\t%d\t\
			%lf\t\
			%*f\t%*d\t%*d\t%*d\t%*f\t%*f\t%*f\t", // these are not assigned to anything
			&temp_atom.type,
			&temp_atom.cart_coord[0], &temp_atom.cart_coord[1], &temp_atom.cart_coord[2],
			&temp_atom.lattice[0], &temp_atom.lattice[1], &temp_atom.lattice[2],
			&temp_atom.bsradius);

		for (j=0;j<MAXIMUM_NUMBER_OF_NEIGHBORS+DISSOLUTION;++j) //when do we pick lattice?
			fscanf(input_file, "%d\t", &temp_atom.transition_indices[j]);

		for (j=0;j<MAXIMUM_NUMBER_OF_NEIGHBORS;++j)
			fscanf(input_file, "%d\t", &temp_atom.neighbor_atom_idxs[j]);
			
		fscanf(input_file, "%d\t%d\t", &temp_atom.next_atom, &temp_atom.previous_atom);

		for (j=0;j<MAXIMUM_NUMBER_OF_COSMETIC_BONDS;++j)
			fscanf(input_file, "%d\t", &tempint);

		fscanf(input_file, "%d\t", &tempint);

		fscanf(input_file, "%lf\t", &(tempdouble[0][0]));

		for (j=0;j<3;++j)
			for (k=0;k<3;++k)
				fscanf(input_file, "%lf\t", &(tempdouble[j][k]));
		
		fscanf(input_file, "%lf\t%lf\t%lf\n", tempdouble[0], tempdouble[1], tempdouble[2]);

		x = temp_atom.lattice[0];
		y = temp_atom.lattice[1];
		z = temp_atom.lattice[2];

		if (atom_at(x, y, z, ss->atom_arr, ss->zone_arr, se) == -1)
		{
			j = add_atom(x,y,z,temp_atom.type, SPECIFIED, ss, se);
		}
	}

	primitive_basis2ucell_params(se->primitive_basis, se->ucell_params);
	// organize(ss->atom_arr, ss->atom_cnt, se->primitive_basis);

	fclose(input_file);
	return true;
}

/*******************************************************************************
*******************************************************************************/

bool process_xyz_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls)
{
	// processes file with .xyz format (number of atoms / comment / type x y z)

	int i;

	char xyz_type[BUFFER_SIZE];
	char* typenames[7]; //can have up to 7 atom types
	int ntypes = 0;
	double xyz_pos[3] = {0.0, 0.0, 0.0};
	double radius;
	int atype;
	char command_string[1024];

	//set_primitive_basis(SC); //is this always true? this should be set somewhere else (beforehand or after?)

	//first line should be the number of atoms
	int nremain; //number of expected atoms

	if (fgets(command_string, BUFFER_SIZE, input_file) == NULL)	// EOF, bad
    {
		fclose(input_file);
		return false;
	}

	nremain = atoi(command_string); //first line of a file is the number of expected atoms

	fgets(command_string, BUFFER_SIZE, input_file); //comment line

	//copy command string here to save comment
	int argsread;

	for (; nremain > 0; --nremain){
		if (fgets(command_string, BUFFER_SIZE, input_file) == NULL)	// EOF
      	{
			fclose(input_file);
			fprintf(temp_log, "ERROR! Ran into EOF, expected %d atoms remaining\n", nremain);
			//organize(atom, atom_cnt); //do I need to call this?
			return false;
		}

		i = ss->atom_cnt;	
				
		create_default_atom(i, ss->atom_arr, se);
		++ss->atom_cnt;
		if (ss->atom_cnt > se->max_atoms)
		{
			fprintf(stderr, "Number of atoms (%lld) is exceeding set maximum (%lld)", ss->atom_cnt, se->max_atoms);
			clean_and_exit(errno);
		}

		if ((argsread = sscanf(command_string, "%s %lf %lf %lf %lf", xyz_type, xyz_pos, xyz_pos+1, xyz_pos+2, &radius)) != 5)
		{
			fprintf(temp_log, "ERROR! Failed to read 4 arguments in .xyz file, only read %d\n", argsread);
			fclose(input_file);
			return false;
        }

		ss->atom_arr[i]->cart_coord[0] = xyz_pos[0];
		ss->atom_arr[i]->cart_coord[1] = xyz_pos[1];
		ss->atom_arr[i]->cart_coord[2] = xyz_pos[2];
		atype = match_atom_type(xyz_type, typenames, &ntypes, temp_log);

		if (atype == -1) //check to see if atom type is successfully added
		{
			for (int i = 0; i < ntypes; ++i){
				free(typenames[i]);
				typenames[i] = NULL;
			}
			//organize(atoms, atom_cnt) //???
			return false;
		}

		ss->atom_arr[i]->type = atype;
		ss->atom_arr[i]->bsradius = radius;

		//vecmul(atom[i]->cart_coord, invert_primitive_basis, atom[i]->lattice); // TODO: need to do this later now!
		
	}

	//TODO: does this now just get turned into the atom name array?
	for (int i = 0; i < ntypes; ++i)
	{
		free(typenames[i]);
		typenames[i] = NULL;
	}
	
	fprintf(temp_log, "Successfully read %lld atoms from .xyz file\n", ss->atom_cnt);
	fclose(input_file);
	//organize(atom, atom_cnt); //???
	return true;
}

/*******************************************************************************
*******************************************************************************/

int match_atom_type(char* type, char* types[], int* num_types, FILE* temp_log) {
	for (int i = 0; i < *num_types; ++i) {
		if (strcmp(type, types[i]) == 0)
			return i+1; //types should start at 1 to be consistent with everything else!
	}
	//must be a new type
	if (*num_types < 7) //change to allow for more discrete atom types?
	{
		//this is a new type and we have space for it
		types[*num_types] = strdup(type);
		(*num_types)++;
		return (*num_types);
	}
	else
	{
		//unrecognized!
		fprintf(temp_log, "ERROR! Did not recognize atom type %s\n", type);
		return -1;
	}
}

/*******************************************************************************
*******************************************************************************/

bool process_kmx_file(FILE* temp_log, FILE* input_file, struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls) {
	int newnat;
	int x,y,z;

	fscanf(input_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&se->primitive_basis[0][0], &se->primitive_basis[0][1], &se->primitive_basis[0][2], 
		&se->primitive_basis[1][0], &se->primitive_basis[1][1], &se->primitive_basis[1][2], 
		&se->primitive_basis[2][0], &se->primitive_basis[2][1], &se->primitive_basis[2][2]); 
			
	fscanf(input_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&se->rmat[0][0], &se->rmat[0][1], &se->rmat[0][2], 
		&se->rmat[1][0], &se->rmat[1][1], &se->rmat[1][2], 
		&se->rmat[2][0], &se->rmat[2][1], &se->rmat[2][2]); 

	fscanf(input_file, "%d %d %d", &se->system_size_x, &se->system_size_y, &se->system_size_z);

	fscanf(input_file, "%d", &newnat);

	fprintf(temp_log, "system size %d %d %d, number of atoms %d\n", se->system_size_x, se->system_size_y, se->system_size_z, newnat);
		
	Atom temp_atom;
	for (int i=0; i<newnat; ++i)
	{
		fscanf(input_file, "%s\t", temp_atom.name);

		fscanf(input_file, 
			"%c\t\
			%lf\t%lf\t%lf\t\
			%d\t%d\t%d\t",
			&temp_atom.type,
			&temp_atom.cart_coord[0], &temp_atom.cart_coord[1], &temp_atom.cart_coord[2],
			&temp_atom.lattice[0], &temp_atom.lattice[1], &temp_atom.lattice[2]); //get rid of lattice coords too?

		for (int j=0; j<se->max_neighbors + se->dissolution; ++j) //when do we pick lattice?
			fscanf(input_file, "%d\t", &temp_atom.transition_indices[j]);

		for (int j=0; j<se->max_neighbors; ++j)
			fscanf(input_file, "%d\t", &temp_atom.neighbor_atom_idxs[j]);
				
		fscanf(input_file, "%d\t%d\t", &temp_atom.next_atom, &temp_atom.previous_atom);

		x = temp_atom.lattice[0];
		y = temp_atom.lattice[1];
		z = temp_atom.lattice[2];

		if (atom_at(x, y, z, ss->atom_arr, ss->zone_arr, se) == -1)
		{
			add_atom(x,y,z,temp_atom.type, SPECIFIED, ss, se);
		}
	}

	primitive_basis2ucell_params(se->primitive_basis, se->ucell_params);
	// organize(ss->atom_arr, ss->atom_cnt, se->primitive_basis);

	fclose(input_file);
		
	return true;
}

/*******************************************************************************
*******************************************************************************/

bool output_log_file(FILE* sim_log_file, int frame_num, double elapsed_stime, double temperature, double overpotential, int atom_cnt, double total_internal_energy)
{
	fprintf(sim_log_file, "![%d]\t", frame_num);
	fprintf(sim_log_file, "time = %lf [s]\ttemperature = %lf [K]\tpotential = %lf [eV]\t", elapsed_stime, temperature, overpotential); // TODO: add iteration number to this
	fprintf(sim_log_file, "atoms = %d\tinternal energy = %lf [eV]\n", atom_cnt, total_internal_energy);
	fflush(sim_log_file);
	return true;
}

bool write_xyz_file(char* xyz_filename, int frame_num, struct SimulationState* ss, struct SimulationEnv* se)
{
	bool is_extended = 1;

	char filename_full[BUFFER_SIZE];
	sprintf(filename_full, "%s_%d.xyz", xyz_filename, frame_num);
	FILE* file = fopen(filename_full, "w+");
	if (file == NULL)
	{
		printf("ERROR! Couldn't open output file %s\n", filename_full);
		fprintf(stderr, "Couldn't open file %s: %s", filename_full, strerror(errno));
		clean_and_exit(errno);
	}

	/* format:
		[number of atoms]
		[comment line - exactly one line]
		[element] [x] [y] [z]
	*/
	// using extended XYZ format
	// (https://docs.ovito.org/reference/file_formats/input/xyz.html#file-formats-input-xyz-extended-format)
	fprintf(file, "%lld\n", ss->atom_cnt); //start with number of atoms
	// TODO: move calculation to initialization fxn
	// TODO: make a simulation_basis or simulation_cell variable
	if (is_extended)
	{
		// 3x3 matrix - rows are cell vectors [preferred]
		int u_range = se->simbox_limits_lat[0][1] - se->simbox_limits_lat[0][0];
		int v_range = se->simbox_limits_lat[1][1] - se->simbox_limits_lat[1][0];
		int w_range = se->simbox_limits_lat[2][1] - se->simbox_limits_lat[2][0];

		// ENHANCE: this would be prettier if se->primitive_basis were transposed'
		// [[u1 u2 u3], [v1 v2 v3], [w1 w2 w3]] vs [[u1 v1 w1], [u2 v2 w2], [u3 v3 w3]]
		int u_lat[3] = {u_range, 0, 0};
		int v_lat[3] = {0, v_range, 0};
		int w_lat[3] = {0, 0, w_range};
		double u_cart[3], v_cart[3], w_cart[3];
		
		lattice2cartesian(u_lat, se->primitive_basis, u_cart);
		lattice2cartesian(v_lat, se->primitive_basis, v_cart);
		lattice2cartesian(w_lat, se->primitive_basis, w_cart);
		
		fprintf(
			file, 
			"Lattice=\"%lf %lf %lf %lf %lf %lf %lf %lf %lf\" ",
			u_cart[0], u_cart[1], u_cart[2],
			v_cart[0], v_cart[1], v_cart[2],
			w_cart[0], w_cart[1], w_cart[2]
		);

		int o_lat[] = {
			se->simbox_limits_lat[0][0],
			se->simbox_limits_lat[1][0],
			se->simbox_limits_lat[2][0]
		};
		double o_cart[3];

		lattice2cartesian(o_lat, se->primitive_basis, o_cart);

		fprintf(
			file, 
			"Origin=\"%lf %lf %lf\" ",
			o_cart[0], o_cart[1], o_cart[2]
		);
		fprintf(file, "pbc=\"T T T\" ");
		fprintf(file, "Properties=id:I:1:species:S:1:pos:R:3 ");
	}
	fprintf(file, "time=%lf temperature=%lf potential=%lf energy=%lf\n", ss->elapsed_stime, ss->temperature, ss->overpotential, ss->total_internal_energy); //need to compute energy here!

	Atom **atoms = ss->atom_arr;
	for (int i = 0; i < ss->atom_cnt; ++i)
	{
		lattice2cartesian(atoms[i]->lattice, se->primitive_basis, atoms[i]->cart_coord);
		fprintf(file, "%d %s %lf %lf %lf %lf\n", i, atoms[i]->name, atoms[i]->cart_coord[0], atoms[i]->cart_coord[1], atoms[i]->cart_coord[2], atoms[i]->bsradius); //name is now element type
	}
	//ball and stick or space filling?
	fclose(file);
	return true;
}
