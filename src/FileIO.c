#include "stdafx.h"
//#include "Mesosim Resources.h"
#include "Defs.h"
#include "Geometry.h"
#include "Prototypes.h"
#include "Vector.h"
#include "Global_Externs.h"
#include "Simulation_Global_Externs.h"
#include "Random.h"

bool get_input_file(char* filename)
{
	char extension[4];
	char* file_ender = strrchr(filename, '.'); //everything after the final '.' in the filename

	if (file_ender == NULL)
	{
		if (sim_log_file == NULL)
			fprintf(temp_log, "ERROR! extension not found in file: %s\n", filename); //for reading arguments
		else
			fprintf(sim_log_file, "ERROR! extension not found in file: %s\n", filename); //should only happen when restarting/checkpointing
	}
	file_ender++;
	extension[3] = '\0';
	strncpy(extension, file_ender,3);

	if (strncmp(extension, "xyz", 3) == 0)
	{
		// open simple x,y,z,type coordinate file
		return process_xyz_file(filename);
	}
	else if (strncmp(extension, "kmc", 3) == 0)
	{
		// open kmc type file
		return process_kmc_file(filename);
	}
	else if (strncmp(extension, "in", 2) == 0)
	{
		//open and process parameter input file
		return process_in_file(filename);
    }
	else if (strncmp(extension, "kmx", 3) == 0)
	{
		//open and process new kmc input type that removes fluff (does this need to happen)
		return process_kmx_file(filename);
    }
	else
	{
		if (sim_log_file == NULL)
			fprintf(temp_log, "ERROR! file extension not recognized: .%s\n", extension); //for reading arguments
		else
			fprintf(sim_log_file, "ERROR! file extension not recognized: .%s\n", extension); //should only happen when restarting/checkpointing
		return false;
    }
}


/*******************************************************************************
*******************************************************************************/

bool process_in_file(char* in_filename) {
	FILE* fileid;
	char parameter_line[200];
	int errnum;
	if ((fileid = fopen(in_filename, "r")) == NULL)
	{
		fprintf(temp_log, "ERROR! Couldn't read .in file %s\n", in_filename);
		return false;
	}
	
	while (fgets(parameter_line, 200, fileid) != NULL) {
		if (strncmp(parameter_line, "restart", 7) == 0) {
			//TODO: restart the simulation from a log file and don't do the rest of the loop
		}
		errnum = parse_input(parameter_line);
		if (errnum != NO_INPUT_ERROR)
		{
			//should write to the temp
			fprintf(temp_log, "ERROR! Had issue reading the following line: \"%s\"\n", parameter_line);
			fclose(fileid);
			return false;
		}
	}
	fclose(fileid);
	return true;
}

/*******************************************************************************
*******************************************************************************/

int parse_input(char* line)
{
	//TODO: use strtok?
	//printf("Trying to parse this line! \"%s\"\n", line);
	char* ptr = line; // line from file
	char cmd[200]; // 
	char params[200]; // parameters parsed from line

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
		if ((argsread = sscanf(params, "%d %d %d", &ssx, &ssy, &ssz)) != 3)
		{
			fprintf(temp_log, "ERROR! Could not correctly read system size parameters %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "temp", 4) == 0) {
		// set the system temperature
		if ((argsread = sscanf(params, "%lf", &temperature)) != 1)
		{
			fprintf(temp_log, "ERROR! Could not correctly read temperature parameter %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
    }
	else if (strncmp(cmd, "seed", 4) == 0) {
		// set the random seed
		if (strncmp(params, "random", 6) == 0) {
			//seed should be based on time
			time_t seedtime;
			time(&seedtime);
			seed = (long int)seedtime;
			fprintf(temp_log, "Using random time seed %ld\n", seed);
		}
		else if (strncmp(params, "default", 7) == 0) {
			seed = DEFAULT_SEED;
			fprintf(temp_log, "Using default time seed %ld\n", seed);
		}
		else if ((argsread = sscanf(params, "%ld", &seed)) != 1) //read in a long int
		{
			fprintf(temp_log, "ERROR! Could not correctly read random seed parameter %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "potential", 9) == 0) {
		//set the potential sweep
		argsread = sscanf(params, "%lf %lf %lf", &initialoverpotential, &overpotential_ramp_rate, &maxoverpotential);
		if (argsread == 1) {
			//constant overpotential
			overpotential_ramp_rate = 0.0;
			maxoverpotential = initialoverpotential;
			fprintf(temp_log, "Using constant potential %lf\n", initialoverpotential);
		}
		else if (argsread != 3) //
		{
			fprintf(temp_log, "ERROR! Could not correctly read potential sweep parameters %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "nne", 3) == 0) {
		// TODO: generalize - refer to spparks?
		//set the nearest neighbor energies - returns number of formats read?
		argsread = sscanf(params, "%lf %lf %lf %lf %lf %lf", nnE, nnE+1 ,nnE+2, nnE+3, nnE+4, nnE+5);
		if (argsread == 1) {
			//something with only A atoms!
			fprintf(temp_log, "Read nearest-neighbor energies for unary system\n");
			for (int i = 1; i <= 5; ++i)
				nnE[i] = 0.;
		}
		if (argsread == 3) {
			//something with only A and B atoms!
			fprintf(temp_log, "Read nearest-neighbor energies for binary system\n");
			nnE[3] = nnE[2];
			nnE[2] = 0.; //ordering for nearest neighbor energies
			nnE[4] = 0.;
			nnE[5] = 0.;
		}
		else if (argsread == 6) {
			//something with A, B, and C atoms!
			fprintf(temp_log, "Read nearest-neighbor energies for ternary system\n");
		}
		else {
			fprintf(temp_log, "ERROR! Read invalid number (%d) of nearest-neighbor energies %s\n", argsread, params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "2nne", 4) == 0) {
		//set the nearest neighbor energies
		argsread = sscanf(params, "%lf %lf %lf %lf %lf %lf", nnnE, nnnE+1 ,nnnE+2, nnnE+3, nnnE+4, nnnE+5);
		if (argsread == 1) {
			//something with only A atoms!
			fprintf(temp_log, "Read 2nd nearest-neighbor energies for unary system\n");
			for (int i = 1; i <= 5; ++i)
				nnnE[i] = 0.;
		}
		if (argsread == 3) {
			//something with only A and B atoms!
			fprintf(temp_log, "Read 2nd nearest-neighbor energies for binary system\n");
			nnnE[3] = nnnE[2];
			nnnE[2] = 0.; //ordering for nearest neighbor energies
			nnnE[4] = 0.;
			nnnE[5] = 0.;
		}
		else if (argsread == 6) {
			//something with A, B, and C atoms!
			fprintf(temp_log, "Read 2nd nearest-neighbor energies for ternary system\n");
		}
		else {
			fprintf(temp_log, "ERROR! Read invalid number (%d) of 2nd nearest-neighbor energies %s\n", argsread, params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "time", 4) == 0) {
		//set time increments for data and max simulation time
		if (strncmp(params, "linear", 6) == 0)
		{
			//read parameters for linear data recording
			argsread = sscanf(params, "%*s %lf %lf %lf", &time_interval_end, &data_time_interval, &run_time);
			analysis_type = REGULAR_TIME_INTERVALS;
		}
		else if (strncmp(params, "log", 3) == 0)
		{
			//read parameters for logarithmic data recording
			argsread = sscanf(params, "%*s %lf %lf %lf", &initial_logtime, &logtime_multiplier, &run_time);
			//initial logtime instead of time_interval_end?
			analysis_type = LOG_TIME_INTERVALS;
		}
		else {
			fprintf(temp_log, "ERROR! Invalid frequency found for data recording %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
		if (argsread != 3) {
			fprintf(temp_log, "ERROR! Could not correctly read time parameters %s\n", params);
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
			lattice_type = FCC;
		else if (strncmp(structtype, "BCC", 3) == 0 || strncmp(structtype, "2", 1) == 0)
			lattice_type = BCC;
		else if (strncmp(structtype, "SC", 2) == 0 || strncmp(structtype, "3", 1) == 0)
			lattice_type = SC;
		/*else if (strncmp(structtype, "DIA", 3) == 0 || strncmp(structtype, "4", 1) == 0)
			lattice_type = DIAMOND;*/ 
		else {
			fprintf(temp_log, "ERROR! Structure type %s not valid\n", structtype);
			return FILE_COMMAND_IGNORED;
		}
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
			simulation_type = SIMULATION_TYPE_FLAT_SHEET;
			if ((argsread = sscanf(params, "%*s %d", &sheet_thickness)) != 1) {
				fprintf(temp_log, "ERROR! Could not correctly read sheet thickness\n");
				return FILE_COMMAND_IGNORED;
			}
		}
		else if (strncmp(params, "cluster", 6) == 0) {
			simulation_type = SIMULATION_TYPE_CLUSTER;
			if ((argsread = sscanf(params, "%*s %d", &cluster_radius)) != 1) {
				fprintf(temp_log, "ERROR! Could not correctly read cluster radius\n");
				return FILE_COMMAND_IGNORED;
			}
		}
		else if (strncmp(params, "file", 4) == 0) {
			simulation_type = SIMULATION_TYPE_FROM_FILE;
			if ((argsread = sscanf(params, "%*s %s", atoms_filename)) != 1) {
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
		char type1[3];
		char type2[3];
		char type3[3];
		argsread = sscanf(params, "%s %s %s", type1, type2, type3);
		switch(argsread) {
			case 1:
				strcpy(atom_names[0], type1);
				break;
			case 2:
				strcpy(atom_names[0], type1);
				strcpy(atom_names[1], type2);
				break;
			case 3:
				strcpy(atom_names[0], type1);
				strcpy(atom_names[1], type2);
				strcpy(atom_names[2], type3);
				break;
			default:
				fprintf(temp_log, "ERROR! Couldn't read any atom type names %s\n", params);
				return FILE_COMMAND_IGNORED;
		}
	}
	else if (strncmp(cmd, "dissolution", 11) == 0) {
		//this determines which atoms dissolve
		char type1[10];
		char type2[10];
		char type3[10];
		int soluble[3];
		argsread = sscanf(params, "%s %s %s", type1, type2, type3);
		if (argsread < 1 || argsread > 3) {
			fprintf(temp_log, "ERROR! Could not correctly read solubilities %s\n", params);
			return FILE_COMMAND_IGNORED;
		}

		soluble[0] = parse_boolean(type1);
		if (soluble[0] == -1) {
			fprintf(temp_log, "ERROR! Could not correctly read solubility %s\n", type1);
			return FILE_COMMAND_IGNORED;
		}
		solubility[0] = soluble[0];

		if (argsread > 1) {
			//only for 2ary and 3ary systems
			soluble[1] = parse_boolean(type2);
			if (soluble[1] == -1) {
				fprintf(temp_log, "ERROR! Could not correctly read solubility %s\n", type2);
				return FILE_COMMAND_IGNORED;	
			}
			solubility[1] = soluble[1];
		}
		if (argsread == 3) {
			soluble[2] = parse_boolean(type3);
			if (soluble[2] == -1) {
				fprintf(temp_log, "ERROR! Could not correctly read solubility %s\n", type3);
				return FILE_COMMAND_IGNORED;
			}
			solubility[2] = soluble[2];
		}
	}
	else if (strncmp(cmd, "composition", 11) == 0) {
		//determines the composition of the atoms: must add up to 1
		double comps[3];
		argsread = sscanf(params, "%lf %lf %lf", comps, comps+1, comps+2);
		if (argsread == 3) {
			if (comps[0] + comps[1] + comps[2] > 1. + 1e-4 || comps[0] + comps[1] + comps[2] < 1. - 1e-4)
			{
				fprintf(temp_log, "ERROR! Compositions don't sum to 1 (%lf + %lf + %lf = %lf)\n", comps[0], comps[1], comps[2], comps[0] + comps[1] + comps[2]);
				return FILE_COMMAND_IGNORED;
			}
			substrate_percent_a = comps[0];
			substrate_percent_b = comps[1];
		}
		else if (argsread == 2) {
			if (comps[0] + comps[1] > 1. + 1e-4 || comps[0] + comps[1] < 1. - 1e-4) {
				fprintf(temp_log, "ERROR! Compositions don't sum to 1 (%lf + %lf = %lf)\n", comps[0], comps[1], comps[0] + comps[1]);
				return FILE_COMMAND_IGNORED;
			}
			substrate_percent_a = comps[0];
			substrate_percent_b = comps[1];
		}
		else if (argsread == 1) {
			//if this isn't the number 1 I will be worried
			if (comps[0] > 1. + 1e-4 || comps[0] < 1. - 1e-4) {
				fprintf(temp_log, "ERROR! Unary system does not have composition 1 (%lf)\n", comps[0]);
				return FILE_COMMAND_IGNORED;
			}
			substrate_percent_a = 1.;
			substrate_percent_b = 0.;
		}
		else {
			fprintf(temp_log, "ERROR! Could not correctly read composition %s\n", params);
			return FILE_COMMAND_IGNORED;
		}
	}
	else
	{
		fprintf(temp_log, "ERROR! keyword %s not recognized\n", cmd);
		return FILE_COMMAND_IGNORED;
	}
	return NO_INPUT_ERROR;
}

int parse_boolean(char *str) {
	if (strncmp(str, "true", 4) == 0 || strncmp(str, "True", 4) == 0 || strncmp(str, "TRUE", 4) == 0 || strncmp(str, "T", 1) == 0 || strncmp(str, "1", 1) == 0)
		return 1;
	else if (strncmp(str, "false", 5) == 0 || strncmp(str, "False", 5) == 0 || strncmp(str, "FALSE", 5) == 0 || strncmp(str, "F", 1) == 0 || strncmp(str, "0", 1) == 0)
		return 0;
	else
		return -1;
}

/*******************************************************************************
*******************************************************************************/

bool process_kmc_file(char *kmc_filename)
	{
		FILE *view_command_file;
		int newnat;

		int i,j,k;
		double x,y,z;

		// [ ]: What are these?
		ssx = DSIMSIZE; //is this always true?????
		ssy = DSIMSIZE;
		ssz = DSIMSIZE;
		zix = TTS;
		ziy = TTS;
		ziz = TTS;							// zones in x, y, z

		if ((view_command_file = fopen(kmc_filename, "r")) == NULL)
			{
				fprintf(temp_log, "ERROR! Couldn't read .kmc file %s\n", kmc_filename);
				return false;
			}

		//general_simulation_initialization(); //happens later

		//read in the lattice and rotation matrices
		fscanf(view_command_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&latmat[0][0], &latmat[0][1], &latmat[0][2], 
			&latmat[1][0], &latmat[1][1], &latmat[1][2], 
			&latmat[2][0], &latmat[2][1], &latmat[2][2]); 
			
		fscanf(view_command_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&rmat[0][0], &rmat[0][1], &rmat[0][2], 
			&rmat[1][0], &rmat[1][1], &rmat[1][2], 
			&rmat[2][0], &rmat[2][1], &rmat[2][2]); 

		fscanf(view_command_file, "%d %d %d", &ssx, &ssy, &ssz);

		fscanf(view_command_file, "%d", &newnat);

		fprintf(temp_log, "system size %d %d %d, number of atoms %d\n", ssx, ssy, ssz, newnat);
		
		int tempint;
		double tempdouble[3][3];

		for (i=0;i<newnat;++i)
		{
			fscanf(view_command_file, "%s\t",
				temp_atom.name);

			fscanf(view_command_file, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%*lf\t%*d\t%*d\t%*d\t%*lf\t%*lf\t%*lf\t",
				&temp_atom.type,
				&temp_atom.coord[0], &temp_atom.coord[1], &temp_atom.coord[2],
				&temp_atom.lattice[0], &temp_atom.lattice[1], &temp_atom.lattice[2],
				&temp_atom.bsradius);

			for (j=0;j<MAXIMUM_NUMBER_OF_NEIGHBORS+DISSOLUTION;++j) //when do we pick lattice?
				fscanf(view_command_file, "%d\t", &temp_atom.position_on_transition_list[j]);

			for (j=0;j<MAXIMUM_NUMBER_OF_NEIGHBORS;++j)
				fscanf(view_command_file, "%d\t", &temp_atom.occupied_neighbor_sites[j]);
				
			fscanf(view_command_file, "%d\t%d\t", &temp_atom.next_atom, &temp_atom.previous_atom);

			for (j=0;j<MAXIMUM_NUMBER_OF_COSMETIC_BONDS;++j)
				fscanf(view_command_file, "%*d\t", &tempint);

			fscanf(view_command_file, "%*d\t", &tempint);
	
			fscanf(view_command_file, "%*lf\t", tempdouble);

			for (j=0;j<3;++j)
				for (k=0;k<3;++k)
					fscanf(view_command_file, "%*lf\t", tempdouble[j][k]);
			
			fscanf(view_command_file, "%*lf\t%*lf\t%*lf\n", tempdouble[0], tempdouble[1], tempdouble[2]);

			x = temp_atom.lattice[0];
			y = temp_atom.lattice[1];
			z = temp_atom.lattice[2];

			if (atom_at(x, y, z) == -1)
			{
				j = add_atom(x,y,z,temp_atom.type, SPECIFIED);
			}
		}

		latmat_to_cell(latmat, cell);
		organize(atom, nat);


		fclose(view_command_file);
		return true;
	}

/*******************************************************************************
*******************************************************************************/

bool process_xyz_file(char *xyz_filename)
{
	// processes file with .xyz format (number of atoms / comment / type x y z)

	int i;
	FILE *view_command_file;

	char xyz_type[200];
	char* typenames[7]; //can have up to 7 atom types
	int ntypes = 0;
	double xyz_pos[3] = {0.0, 0.0, 0.0};
	double radius;
	int atype;

	if ((view_command_file = fopen(xyz_filename, "r")) == NULL)
	{
		fprintf(temp_log, "ERROR! Couldn't read .xyz file %s\n", xyz_filename);
		return false;
	}

	//set_latmat(SC); //is this always true? this should be set somewhere else (beforehand or after?)

	//first line should be the number of atoms
	int nremain; //number of expected atoms

	if (fgets(command_string, 200, view_command_file) == NULL)	// EOF, bad
    {
		fclose(view_command_file);
		return false;
	}

	nremain = atoi(command_string); //first line of a file is the number of expected atoms

	fgets(command_string, 200, view_command_file); //comment line

	//copy command string here to save comment
	int argsread;

	for (; nremain > 0; --nremain){
		if (fgets(command_string, 200, view_command_file) == NULL)	// EOF
      	{
			fclose(view_command_file);
			fprintf(temp_log, "ERROR! Ran into EOF for %s, expected %d atoms remaining\n", xyz_filename, nremain);
			//organize(atom, nat); //do I need to call this?
			return false;
		}

		i=nat;	
				
		create_default_atom(i);

		if ((argsread = sscanf(command_string, "%s %lf %lf %lf %lf", xyz_type, xyz_pos, xyz_pos+1, xyz_pos+2, &radius)) != 5)
		{
			fprintf(temp_log, "ERROR! Failed to read 4 arguments in .xyz file, only read %d\n", argsread);
			fclose(view_command_file);
			return false;
        }

		atom[i]->coord[0] = xyz_pos[0];
		atom[i]->coord[1] = xyz_pos[1];
		atom[i]->coord[2] = xyz_pos[2];
		atype = match_atom_type(xyz_type, typenames, &ntypes);

		if (atype == -1) //check to see if atom type is successfully added
		{
			for (int i = 0; i < ntypes; ++i)
				free(typenames[i]);
			//organize(atoms, nat) //???
			return false;
		}

		atom[i]->type = atype;
		atom[i]->bsradius = radius;

		//vecmul(atom[i]->coord, ilatmat, atom[i]->lattice); // TODO: need to do this later now!

        /*switch(atom[i]->type)
        {
			case 0:  				// white
				atom[i]->color[0] = 1.;
   				atom[i]->color[1] = 1.;
      		   	atom[i]->color[2] = 1.;
				break;

			case 1:					// red	
		   	    atom[i]->color[0] = 1.;
	   			atom[i]->color[1] = 0.;
		  		atom[i]->color[2] = 0.;
               	break;

			case 2:					// green
			   	atom[i]->color[0] = 0.;
	   			atom[i]->color[1] = 1.;
      			atom[i]->color[2] = 0.;
				break;
			
			case 3:					// blue
				atom[i]->color[0] = 0.;
   				atom[i]->color[1] = 0.;
      			atom[i]->color[2] = 1.;
				break;

			case 4:					// cyan
				atom[i]->color[0] = 0.;
			   	atom[i]->color[1] = 1.;
      			atom[i]->color[2] = 1.;
		        break;
        
			case 5:					// gray
				atom[i]->color[0] = .5;
				atom[i]->color[1] = .5;
				atom[i]->color[2] = .5;
				break;

			case 6:					// yellow
				atom[i]->color[0] = 1.;
				atom[i]->color[1] = 1.;
				atom[i]->color[2] = 0.;
				break;
	    }*/

		++nat;
		
	}

	//TODO: does this now just get turned into the atom name array?
	for (int i = 0; i < ntypes; ++i)
		free(typenames[i]);
	
	fprintf(temp_log, "Successfully read %d atoms from .xyz file %s\n", nat, xyz_filename);
	fclose(view_command_file);
	//organize(atom, nat); //???
	return true;
}

/*******************************************************************************
*******************************************************************************/

int match_atom_type(char* type, char* types[], int* num_types) {
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

bool process_kmx_file(char* kmx_filename) {
	FILE *view_command_file;

	int newnat;
	int i,j,k;
	double x,y,z;


	if ((view_command_file = fopen(kmx_filename, "r")) == NULL)
	{
		fprintf(temp_log, "ERROR! Couldn't read .kmx file %s\n", kmx_filename);
		return false;
	}

	fscanf(view_command_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&latmat[0][0], &latmat[0][1], &latmat[0][2], 
		&latmat[1][0], &latmat[1][1], &latmat[1][2], 
		&latmat[2][0], &latmat[2][1], &latmat[2][2]); 
			
	fscanf(view_command_file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&rmat[0][0], &rmat[0][1], &rmat[0][2], 
		&rmat[1][0], &rmat[1][1], &rmat[1][2], 
		&rmat[2][0], &rmat[2][1], &rmat[2][2]); 

	fscanf(view_command_file, "%d %d %d", &ssx, &ssy, &ssz);

	fscanf(view_command_file, "%d", &newnat);

	fprintf(temp_log, "system size %d %d %d, number of atoms %d\n", ssx, ssy, ssz, newnat);
		

	for (i=0;i<newnat;++i)
	{
		fscanf(view_command_file, "%s\t",
			temp_atom.name);

		fscanf(view_command_file, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",
			&temp_atom.type,
			&temp_atom.coord[0], &temp_atom.coord[1], &temp_atom.coord[2],
			&temp_atom.lattice[0], &temp_atom.lattice[1], &temp_atom.lattice[2]); //get rid of lattice coords too?

		for (j=0;j<MAXIMUM_NUMBER_OF_NEIGHBORS+DISSOLUTION;++j) //when do we pick lattice?
			fscanf(view_command_file, "%d\t", &temp_atom.position_on_transition_list[j]);

		for (j=0;j<MAXIMUM_NUMBER_OF_NEIGHBORS;++j)
			fscanf(view_command_file, "%d\t", &temp_atom.occupied_neighbor_sites[j]);
				
		fscanf(view_command_file, "%d\t%d\t", &temp_atom.next_atom, &temp_atom.previous_atom);

		x = temp_atom.lattice[0];
		y = temp_atom.lattice[1];
		z = temp_atom.lattice[2];

		if (atom_at(x, y, z) == -1)
		{
			j = add_atom(x,y,z,temp_atom.type, SPECIFIED);
		}
	}

	latmat_to_cell(latmat, cell);
	organize(atom, nat);


	fclose(view_command_file);
		
	return true;
}

/*******************************************************************************
*******************************************************************************/

bool output_log_file(int frame_num)
{
	fprintf(sim_log_file, "![%d]\t", frame_num);
	fprintf(sim_log_file, "time = %lf [s]\ttemperature = %lf [K]\tpotential = %lf [eV]\t", elapsed_time, temperature, overpotential);
	fprintf(sim_log_file, "atoms = %d\tinternal energy = %lf [eV]\n", nat, total_internal_energy);
	return true;
}

bool write_xyz_file(char* xyz_filename, int frame_num)
{
	FILE* fileid;
	char filename_full[200];
	sprintf(filename_full, "%s_%d.xyz", xyz_filename, frame_num);
	if ((fileid = fopen(filename_full, "w+")) == NULL)
	{
		printf("ERROR! Couldn't open output file %s\n", filename_full);
		return false;
	}
	fprintf(fileid, "%d\n", nat); //start with number of atoms
	fprintf(fileid, "time = %lf, temperature = %lf, potential = %lf, energy = %lf\n", elapsed_time, temperature, overpotential, total_internal_energy); //need to compute energy here!

	for (int i = 0; i < nat; ++i)
		fprintf(fileid, "%s %lf %lf %lf %lf\n", atom[i]->name, atom[i]->coord[0], atom[i]->coord[1], atom[i]->coord[2], atom[i]->bsradius); //name is now element type
	//ball and stick or space filling?
	fclose(fileid);
	return true;
}