#include "stdafx.h"
#include "Defs.h"
#include "Geometry.h"
#include "Vector.h"
#include "Random.h"
#include "Simulation.h"
#include "Simulation_Aux.h"
#include "Atoms.h"
#include "FileIO.h"

// Simulation
int rate_skip;
double xr, yr, zr;
int rwx, rwy, rwz;
double rw[3], rrp[3];
double ta1, ta2;
int adatom_before;
char g_atom_names[3][3]={"1", "2", "3"};
double g_default_color[3] = {0., 0., 0.};
int g_simulation_type = SIMULATION_TYPE_UNDEFINED;
int g_atom_cnt = 0;
Atom g_temp_atom;
// bool ss->simulation_should_kill_itself;
Atom* g_atom_arr[]; // array containing all atoms in the simulation
// double ss->elapsed_stime = 0; 
bool g_evaporation_flag = true;
char g_coordinate_log_prefix[256] = "default_simulation_analysis.dat";
// unsigned long ss->final_iteration = 1e9; // TODO: move to input file
// double ss->run_stime = 1.e8; //default simulation runtime (in seconds) // TODO: move to input file
double g_log_interval = 0.1; // logging rate in simulation time
double g_next_log_checkpoint; // simulation time when next log will be output
double* g_log_list = NULL; // pointer to list of points at which to log data
int g_log_list_len = 0; // length of g_log_list (number of times to log data)
// int ss->sim_end_type = 0;
int g_analysis_type = REGULAR_TIME_INTERVALS;
int g_rate_cnt; // number of rates in g_rate_arr list (filled indices)
int g_transition_cnt; // size of filled portion of transition list
// double ss->frequency_sum;
// double ss->overpotential = 0.0;
double g_nnE[6] = {
	DEFAULT_BOND_ENERGY_AA,
	DEFAULT_BOND_ENERGY_AB,
	DEFAULT_BOND_ENERGY_AC,
	DEFAULT_BOND_ENERGY_BB,
	DEFAULT_BOND_ENERGY_BC,
	DEFAULT_BOND_ENERGY_CC};
double g_nnnE[6] = {0., 0., 0., 0., 0., 0.};
bool g_solubility[3] = {false, false, false}; // whether ABC-type atoms can be dissolved/evaporated - all elements cannot dissolve by default
// double ss->temperature = DEFAULT_TEMPERATURE;
int g_dissolution = DISSOLUTION;
int g_final_config_neighbor_cnt;
int g_intial_config_neighbor_cnt;
// ENHANCE: malloc?
Rate g_rate_arr[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS]; // array containing all the unique rate constants and count of atoms that have that k and indices in g_transition_arr
// Rates default initalized to all zeros values?
Transition *g_transition_arr[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS]; // array of the possible atom transitions (atom index + jump offset index)
//Transition g_transition_arr[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
// contains all the same information in g_atom_arr[i].transition_indices
// Trans_Prob g_transition_probability;
int g_lastxt, g_lastyt, g_lastzt;
double g_sum_of_rate_populations;
double g_current_probability;
bool checkpoint_reached = false;


// Simulation_Aux
// double ss->total_internal_energy = 0;
// [ ]: what are these?
int g_zixshift, g_ziyshift, g_zizshift; // bit shifts for finding zones from coordinates
int g_ssxshift, g_ssyshift, g_sszshift; // used with zi*shift
int g_zsh, g_ysh, g_xsh;	// total bit shifts, zi*shift - ss*shift
// [ ]: what are the units for this? how does it relate to atomic spacing?
int g_ssx = DSIMSIZE, g_ssy = DSIMSIZE, g_ssz = DSIMSIZE;	// system size x, y, z in lattice coordinates // TODO: lattice coordinates but not along lattice vectors? // XXX: overwritten by input file
double g_ssr;
int g_zix = TTS, g_ziy = TTS, g_ziz = TTS;
int g_lattice_type = FCC;
int g_max_neighbors = 12; // [ ]: this should be dependent on the crystal structure
int g_sheet_thickness = -1;
int g_cluster_radius = -1;
char g_atoms_filename[256] = "";
// Zone g_zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z];
double g_initial_overpotential = DEFAULT_OVERPOTENTIAL;
double g_overpotential_ramp_rate = 0.0;
double g_max_overpotential = DEFAULT_OVERPOTENTIAL;
double g_substrate_percent_a = DEFAULT_COMPOSITION_A;
double g_substrate_percent_b = DEFAULT_COMPOSITION_B;
// int ss->total_volume_dissolved;
double g_normal_x, g_normal_y, g_normal_z;
double lhs[6];
double normal_lat[6][3];
int translation_vector[6][3];

// ENHANCE: pass struct with all simulation parameters as argument
unsigned long perform_simulation(struct SimulationState *ss) //potentially FILE* as arguments
{
	//printf("Simulation starting\n");
	long int i;
	int j, k;

	double nt, ot;

	double transition_type_probability;

	//double vap;

	int which_one;
	int atom_number, jump_vector;
	int natn; //changed from nan bc keyword
	int atype;

	int moved_flag = true;

	int framenum = 0; // [ ]: what is this? it isn't iteration count

	ss->elapsed_stime = 0.0;
	//writes data time intervals and run time in original code

	//simulation_is_going = true;

	//print system and zone sizes to the file?

	/*if (g_num_sims > 0) //this probably will not happen
	{
		do_initialize_simulation(g_simulation_type);
		ss->elapsed_stime = 0.0;
	}*/ 

	// initialize simulation kinetics, and draw a picture

	//printf("transition time!\n");
	for (i=0;i<g_atom_cnt;++i)	
		refresh_transitions(i, ss);			// resets all kinetic paramters

	organize(g_atom_arr, g_atom_cnt); //replacement for copy_xyz_to_coord but might not be necessary

	//printf("transitioned and organized\n");

	ss->total_volume_dissolved = 0; //do i care

	// initial state
	calculate_internal_energy(g_atom_cnt, ss);
	//printf("energy calculated\n");
	output_log_file(g_sim_log_file, framenum, ss);
	write_xyz_file(g_coordinate_log_prefix, framenum, ss);
	//printf("files written\n");
	++framenum;

	unsigned long int iter = 1; // iteration count

	ot = 0.0; //needs to happen outside of loop
	//printf("about to start the loop\n");
	
	bool simulation_end = false;

	while (!simulation_end)
	{
		if (iter % 100 == 0)
			printf("iteration %ld, time %lf\n", iter, ss->elapsed_stime);
		if (ss->simulation_should_kill_itself) // abort simulation (only happens if atoms overlap)
		{
			//find_average_curvature(); //no longer valid

			calculate_internal_energy(g_atom_cnt, ss);
			output_log_file(g_sim_log_file, framenum, ss);
			write_xyz_file(g_coordinate_log_prefix, framenum, ss);

			ss->simulation_should_kill_itself = false;
			organize(g_atom_arr, g_atom_cnt); //replace the copy with draw, maybe not needed

			//simulation_is_going = false; //not needed
			return 1; //return 1 b/c error?
		}


		for (j=0;j<g_atom_cnt;++j)	
			refresh_transitions(j, ss);			// resets all kinetic paramters
		 //does this need to happen?

		// increment the elapsed time

		compute_transition_array(ss);
		ss->elapsed_stime -= log(drandj(&rand_seed)) / ss->frequency_sum;

		if (ss->sim_end_type == SIM_END_BY_STIME && ss->elapsed_stime >= ss->run_stime) // simulation has gone past time
			break; //get outta here before I make a new transition

		// pick the type of transition to occur
		transition_type_probability = drandj(&rand_seed);
		rate_skip = g_rate_cnt/2;
		j = rate_skip;

		//change structure to see if the move gets made or not

		// binary search to select transition
		moved_flag = false; //used to track diffusion/evaporation vs deposition
		while (moved_flag == false) {
			if ((transition_type_probability >= ss->transition_probability.lbound[j])
						&& (transition_type_probability < ss->transition_probability.ubound[j]))
			{
				k = ss->transition_probability.listnum[j];

				// pick the lucky atom
				// [ ]: third random number?
				// picks a type of transition and then which atom that has that transition will it act on?
				which_one = g_rate_arr[k].transition_start_idx + (int)(drandj(&rand_seed)*(double)g_rate_arr[k].transition_count);

				// which_one gives the location of the g_transition_arr, which gives
				// the info about the specific atom
				atom_number = g_transition_arr[which_one]->atom_idx;
				jump_vector = g_transition_arr[which_one]->offset_idx;

				// if jump_vector == g_max_neighbors then the atom is going to evaporate
				moved_flag = true;

				adatom_before = 0;

				if (jump_vector != g_max_neighbors) //diffusion
				{
					//printf("the transition is diffusion of atom %d, jumping from %lf %lf %lf ", atom_number, atom[atom_number]->lattice[0], atom[atom_number]->lattice[1], atom[atom_number]->lattice[2]);
					// coordinates atom is jumping to
					g_lastxt = g_atom_arr[atom_number]->lattice[0] + jump_offset[jump_vector].dx;
					g_lastyt = g_atom_arr[atom_number]->lattice[1] + jump_offset[jump_vector].dy;
					g_lastzt = g_atom_arr[atom_number]->lattice[2] + jump_offset[jump_vector].dz;

					adjust_pbc(&g_lastxt, &g_lastyt, &g_lastzt);

					atype = g_atom_arr[atom_number]->type;

					// moves atom?
					remove_atom(atom_number, ss);
					natn = add_atom(g_lastxt, g_lastyt, g_lastzt, atype, NORMAL, ss);
					//printf("and jumping to %lf %lf %lf\n", g_lastxt, g_lastyt, g_lastzt);
				}
				else				// dissolution
				{
					if (g_solubility[g_atom_arr[atom_number]->type - 1] == true)	// atoms are dissolved based on input specs!
					{
						++ss->total_volume_dissolved;
						remove_atom(atom_number, ss);	// evaporate the atom
						//printf("the transition was dissolution of atom %d\n", atom_number);
					}
				}
			}
			else if (transition_type_probability < ss->transition_probability.lbound[j])
			{
				//search to the left
				rate_skip = rate_skip/2;
				if (rate_skip == 0) rate_skip = 1;
				j -= rate_skip;
			}
			else if (transition_type_probability >= ss->transition_probability.ubound[j])
			{
				//search to the right
				rate_skip = rate_skip/2;
				if (rate_skip == 0) rate_skip = 1;
				j += rate_skip;
				if (j == g_rate_cnt) //no more options!
					break;
			}
		}

		if (moved_flag == false) //only happens iff jump_vector == g_max_neighbors
		{
			printf("for some reason I didn't transition\n");
			//deposition is vestigial and we don't want it!
			//what happens if we get to this point then?
			/*if (deposition_type == DEPOSITION_TYPE_RAINFALL)
			{
				get_vapor_deposition_site(&dep_x, &dep_y, &dep_z);

				vap = drandj(&rand_seed);
				if (vap < deposition_rate_of_a/(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c))
					atype = 1; 
				else if (vap < (deposition_rate_of_a + deposition_rate_of_b)/(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c))
					atype = 2;
				else atype = 3;

				natn = add_atom(dep_x, dep_y, dep_z, atype, NORMAL);
			}
			else if (deposition_type == DEPOSITION_TYPE_RANDOM_WALKER)
			{
				vap = drandj(&rand_seed);
				if (vap < deposition_rate_of_a/(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c))
					atype = 1; 
				else if (vap < (deposition_rate_of_a + deposition_rate_of_b)/(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c))
					atype = 2;
				else atype = 3;

				// start random walker from sphere touching edge of system.
				// it will then diffuse according to normal diffusion physics (zero bonds) as given by calcrate

				//newt:
				do
				{
					ta1 = 2.*PI*drandj(&rand_seed);
					ta2 = 2.*PI*drandj(&rand_seed);

					rrp[0] = g_ssr*cos(ta2)*cos(ta1);
					rrp[1] = g_ssr*cos(ta2)*sin(ta1);
					rrp[2] = g_ssr*sin(ta2);

					// now invert xr, yr, zr into lattice vectors;

					vecmul(rrp, invert_primitive_basis, rw);

					rwx = (int)rw[0] + g_ssx/2;
					rwy = (int)rw[1] + g_ssy/2;
					rwz = (int)rw[2] + g_ssz/2;

				} while (atom_at(rwx, rwy, rwz, ss) >= 0);

				//if (atom_at(rwx, rwy, rwz, ss) >= 0) goto newt; (made redundant with do while)

				natn = add_atom(rwx, rwy, rwz, atype, NORMAL);
			}*/
		}
		 
		// after iteration, log if necessary
		// TODO: implement the checkpoint lists
		if ((g_analysis_type == REGULAR_TIME_INTERVALS) || (g_analysis_type == LN_TIME_INTERVALS))
			checkpoint_reached = (ss->elapsed_stime >= g_next_log_checkpoint);
		else if (g_analysis_type == ITERATION_INTERVALS)
			checkpoint_reached = (iter >= g_next_log_checkpoint);

		if (checkpoint_reached)
		{
			organize(g_atom_arr, g_atom_cnt); //replaced but do i really need it
			
			//record the elapsed time in a file here
			printf("writing file %d: elapsed_stime = %lf\n", framenum, ss->elapsed_stime);
			
			calculate_internal_energy(g_atom_cnt, ss);
			output_log_file(g_sim_log_file, framenum, ss);
			write_xyz_file(g_coordinate_log_prefix, framenum, ss);
			
			if (g_analysis_type == REGULAR_TIME_INTERVALS)
			{
				// bring g_next_log_checkpoint up to and one step beyond ss->elapsed_stime
				while (g_next_log_checkpoint <= ss->elapsed_stime)
					g_next_log_checkpoint += g_log_interval;

			}
			else if (g_analysis_type == LN_TIME_INTERVALS)
			{
				while (g_next_log_checkpoint <= ss->elapsed_stime)
					g_next_log_checkpoint *= g_log_interval;
			}
			else if (g_analysis_type == ITERATION_INTERVALS)
			{
				g_next_log_checkpoint += g_log_interval;
			}

			if (g_overpotential_ramp_rate != 0.0)
			{
				nt = ss->elapsed_stime;
				ss->overpotential += (nt-ot)*g_overpotential_ramp_rate;
				ot = ss->elapsed_stime;
				for (j=0;j<g_atom_cnt;++j)	
					refresh_transitions(j, ss);			// resets all kinetic paramters

			}
			
			++framenum;
		}

		++iter; //sanity check to avoid ending in an infinite cycle // [ ]: what?

		// check if simulation is over
		if (ss->sim_end_type == SIM_END_BY_STIME){
			simulation_end = (ss->elapsed_stime >= ss->run_stime);
		}
		else if (ss->sim_end_type == SIM_END_BY_ITERATIONS) {
			simulation_end = (iter >= ss->final_iteration);
		}
	}

	if (iter == ss->final_iteration)
		fprintf(g_sim_log_file, "reached final iteration and terminated\n");
	
	//write ss->elapsed_stime to mark finish
	//TODO: finish IO
	calculate_internal_energy(g_atom_cnt, ss);
	output_log_file(g_sim_log_file, framenum, ss);
	write_xyz_file(g_coordinate_log_prefix, framenum, ss);

	printf("Finished simulation\n"); //move this to the log file
		
	//simulation_is_going = false;
	return 0;
}

/******************************************************************************/
/******************************************************************************/
// updates transition_probabilty (weighted rate list, used to choose event)
void compute_transition_array(struct SimulationState *ss)
{
	int i;
	int x;

	int total_lists = 0;

	ss->frequency_sum = 0.0;
	g_sum_of_rate_populations = 0.0;

	for (i=0;i<g_rate_cnt;++i)
	{
		if (g_rate_arr[i].transition_count != 0)
		{
			g_rate_arr[i].frequency = g_rate_arr[i].k*(double)g_rate_arr[i].transition_count;
			ss->frequency_sum += g_rate_arr[i].frequency;
			g_sum_of_rate_populations += g_rate_arr[i].transition_count;

			ss->transition_probability.listnum[total_lists] = i;
			++total_lists;
		}
	}

	/*if (deposition_type == DEPOSITION_TYPE_RAINFALL)
		ss->frequency_sum += (double)g_ssx*(double)g_ssy*(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c);
	else if (deposition_type == DEPOSITION_TYPE_RANDOM_WALKER)
		ss->frequency_sum += (4*PI*g_ssr*g_ssr)*(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c);*/ //deposition is vestigial

	// now compute bounds for jump probabilities

	g_current_probability = 0.0;

	for (i=0;i<total_lists;++i)
	{
		ss->transition_probability.lbound[i] = g_current_probability;
		x = ss->transition_probability.listnum[i];
		g_current_probability += g_rate_arr[x].frequency/ss->frequency_sum;
		ss->transition_probability.ubound[i] = g_current_probability;
	}

	return;
}

/******************************************************************************/
/******************************************************************************/
// updates [g_atom_arr[atom_idx], g_transition_arr[i], g_rate_arr[i].transition_start_idx], initializes g_rate_arr
int refresh_transitions(int atom_idx, struct SimulationState *ss) // atom_idx = index on atom list
{
	int i, j;
	int rate_idx; // position of rate rate in rate list g_rate_arr[]
	double rate; // rate constant from bond-breaking model
	int next_x, next_y, next_z;

	int atom_rates_cnt;	// this is returned as the number of transitions (jump_offsets) this atom can undergo, excluding evaporation
	// 
	int start_config[MAXIMUM_NUMBER_OF_NEIGHBORS]; // -1 if empty, type if filled
	int end_config[MAXIMUM_NUMBER_OF_NEIGHBORS];

	//bool ok_to_evaporate; //doesn't get used // XXX: commented code

	//printf("removing\n"); // XXX: commented print
	// first, remove all mention of this atom from transition list
	for (i=0; i<g_max_neighbors + g_dissolution; ++i) // extra 1 for evaporation
	{
		if (g_atom_arr[atom_idx]->transition_indices[i] != -1) // something can happen in the "i" direction
			take_off_transition_list(atom_idx, i);
	}

	//printf("cycling\n"); // XXX: commented print
	// cycle through neighbor coordinates, and check if there is an atom there
	for (i=0;i<g_max_neighbors;++i)
	{
   		next_x = g_atom_arr[atom_idx]->lattice[0] + jump_offset[i].dx;
		next_y = g_atom_arr[atom_idx]->lattice[1] + jump_offset[i].dy;
        next_z = g_atom_arr[atom_idx]->lattice[2] + jump_offset[i].dz;

		adjust_pbc(&next_x, &next_y, &next_z);

        j = atom_at(next_x, next_y, next_z, ss);

		if (j >= 0) g_atom_arr[atom_idx]->occupied_neighbor_sites[i] = j;
		// if no atom atom_idx that position but occ_neighbor array says there is, fix it
		if ((j == -1)&&(g_atom_arr[atom_idx]->occupied_neighbor_sites[i] >= 0))
			g_atom_arr[atom_idx]->occupied_neighbor_sites[i] = -1;
	}

	// cycle through the neighbor sites.  if there's an empty one, calculate the transition rate to it

	//printf("getting config\n"); // XXX: commented prints
	atom_rates_cnt = 0;
	//ok_to_evaporate = true; //doesn't get used

	g_intial_config_neighbor_cnt = get_initial_configuration2(atom_idx, 0, start_config);		// k is number of near neighbors

	if (g_intial_config_neighbor_cnt == g_max_neighbors) // skip calculating a rate of a fully coordinated atom
		return atom_rates_cnt;			

	//printf("calculating surf diffusion rate\n");
	for (i=0; i<g_max_neighbors; ++i)
	{
		if (g_atom_arr[atom_idx]->occupied_neighbor_sites[i] == -1)
		{ // if unoccupied, consider transition
			//printf("final config\n");
			g_final_config_neighbor_cnt = get_final_configuration2(atom_idx, i, end_config, ss);
			//printf("surf diffusion\n");
			calculate_surf_diffusion_rate(		start_config, // ENHANCE: make this look prettier
												end_config,
												g_max_neighbors,
												g_atom_arr[atom_idx]->type,
												g_nnE,
												ss->temperature,
												ss->overpotential,
												&rate);
					
			++atom_rates_cnt;
			//printf("adding step\n");
			if ((rate_idx = is_on_transition_list(rate)) != -1)
			{ // if rate constant already in rate list, use that rate list index
				add_to_transition_list(rate_idx, atom_idx, i);
			}
			else
			{
				// the transition rate with rate constant rate turned out to be a new one
				// add to rate list/array
				rate_idx = create_new_transition(rate);
				add_to_transition_list(rate_idx, atom_idx, i);
			}
		}
	}	

	//printf("calcualte evap\n");
	calculate_evaporation_rate(	start_config,
								g_max_neighbors,
								g_atom_arr[atom_idx]->type,
								g_nnE,
								ss->temperature,
								ss->overpotential,
								&rate);								

	//replace with the lines below because it's more obvious
	/*if ((rate_idx = is_on_transition_list(rate)) != -1)
	{
		add_to_transition_list(rate_idx, atom_idx, g_max_neighbors);
	}
	else  // the transition rate to that spot turned out to be a new one.
	{		
		rate_idx = create_new_transition(rate);
		add_to_transition_list(rate_idx, atom_idx, g_max_neighbors);
	}*/

	//printf("end stuff\n");
	rate_idx = is_on_transition_list(rate);

	if (rate_idx == -1)
		rate_idx = create_new_transition(rate); //the transition rate to the spot is a new one!

	add_to_transition_list(rate_idx, atom_idx, g_max_neighbors); // evaporation is considered to be last in jump_offset (not really in array but uses that index number)
	//printf("gonna return %d\n", atom_rates_cnt);
	return atom_rates_cnt;						// gives number of current transitions for that atom
}

/******************************************************************************/
/******************************************************************************/

//is there a better way of doing this?
// checks if the rate constant is already in the rate list
int is_on_transition_list(double rate) 
{
	int i;

	for (i=0;i<g_rate_cnt;++i)
		if (rate == g_rate_arr[i].k) return i;

	return -1;
}

/******************************************************************************/
/******************************************************************************/
// create new Rate struct in rate array
// updates rate_array[g_rate_cnt], g_rate_cnt
int create_new_transition(double rate)
{
	g_rate_arr[g_rate_cnt].k = rate;
	g_rate_arr[g_rate_cnt].transition_start_idx = g_transition_cnt;
	g_rate_arr[g_rate_cnt].transition_count = 0;

	++g_rate_cnt;

	return (g_rate_cnt-1);
}

/******************************************************************************/
/******************************************************************************/

// add to g_rate_arr[rate_idx] the atom atom_idx going in direction offset_idx
// updates g_transition_arr, g_rate_arr[rate_idx].transition_count, g_atom_arr[atom_idx]->transition_indices[offset_idx]
void add_to_transition_list(int rate_idx, int atom_idx, int offset_idx) // g_rate_arr index, g_atom_arr index, jump_offset index
{
	int i; // loop variable
	int n;
	int initial_transition_index, final_transition_index; // initial and final g_transition_arr index

	// make room for the new arrival

	g_transition_arr[g_transition_cnt] = (Transition *)malloc(sizeof(Transition));	// adds entry to the end of the list

	// what is this
	//		final_transition_index = g_rate_arr[g_rate_cnt-1].transition_start_idx + g_rate_arr[g_rate_cnt-1].number;	
	//		g_transition_arr[final_transition_index] = (Transition *)malloc(sizeof(Transition));

	for (i = g_rate_cnt-1;i>rate_idx;--i)
		{ // [ ]: what does this do? is this the same as in remove_transition?
			initial_transition_index = g_rate_arr[i].transition_start_idx;
			final_transition_index = initial_transition_index + g_rate_arr[i].transition_count;

			g_transition_arr[final_transition_index]->atom_idx = g_transition_arr[initial_transition_index]->atom_idx;
			g_transition_arr[final_transition_index]->offset_idx = g_transition_arr[initial_transition_index]->offset_idx;

			if (initial_transition_index != final_transition_index)
				g_atom_arr[g_transition_arr[final_transition_index]->atom_idx]->transition_indices[g_transition_arr[final_transition_index]->offset_idx] = final_transition_index;

			++g_rate_arr[i].transition_start_idx;
		}

	// add new arrival

	n = g_rate_arr[rate_idx].transition_start_idx + g_rate_arr[rate_idx].transition_count;

	++g_rate_arr[rate_idx].transition_count;

	g_transition_arr[n]->atom_idx = atom_idx;
	g_transition_arr[n]->offset_idx = offset_idx;

	g_atom_arr[atom_idx]->transition_indices[offset_idx] = n;

	++g_transition_cnt;

	return;
}

/******************************************************************************/
/******************************************************************************/
// updates g_atom_arr[atom_idx], g_transition_arr[i], g_rate_arr[i].transition_start_idx
void take_off_transition_list(int atom_idx, int offset_idx)	// removes atom jumping in the jump_offset[offset_idx] direction
{
	int i;
	int rate_idx, transition_idx, transition_end_idx;

	// find out what Rate in g_rate_arr this is

	transition_idx = g_atom_arr[atom_idx]->transition_indices[offset_idx]; // old position on transition list, to be removed

	for (i=0;i<g_rate_cnt;++i)
		if (transition_idx < (g_rate_arr[i].transition_start_idx + g_rate_arr[i].transition_count))
		{
			rate_idx = i;
			break;
		}

	// remind atom it can no longer jump

	g_atom_arr[atom_idx]->transition_indices[offset_idx] = -1;
	// TODO: do these later, after it's been removed and g_transition_arr has been rearranged
	--g_transition_cnt;

	// rate_idx points to the current rate list it's on.  decrement the number of atoms in that list
	// and clean up.  If the list is empty, remove it.

	--g_rate_arr[rate_idx].transition_count;
	// [ ]: wtf is going on here 
	if (g_rate_arr[rate_idx].transition_count == 0) // if list is empty
	{
		for (i = rate_idx + 1; i < g_rate_cnt; ++i)
		{
			--g_rate_arr[i].transition_start_idx; // move g_rate_arr offsets of larger indicies down one
			// make transition at new start index (which is of different Rate than old start index) have same atom and offset as new end index (which is of same Rate as old end index)
			transition_idx = g_rate_arr[i].transition_start_idx;
			transition_end_idx = g_rate_arr[i].transition_start_idx+g_rate_arr[i].transition_count;		// count is always at least 1

			g_transition_arr[transition_idx]->atom_idx = g_transition_arr[transition_end_idx]->atom_idx;
			g_transition_arr[transition_idx]->offset_idx = g_transition_arr[transition_end_idx]->offset_idx;
			// update the transition index in the corresponding atom in g_atom_arr to have the new (lower) transition index 
			g_atom_arr[g_transition_arr[transition_idx]->atom_idx]->transition_indices[g_transition_arr[transition_idx]->offset_idx] = transition_idx;
		}

		free(g_transition_arr[g_transition_cnt]);			// free up the very last member of the last g_transition_arr

		for (i=rate_idx+1;i<g_rate_cnt;++i)
		{
			g_rate_arr[i-1].transition_start_idx = g_rate_arr[i].transition_start_idx;
			g_rate_arr[i-1].transition_count = g_rate_arr[i].transition_count;
			g_rate_arr[i-1].k = g_rate_arr[i].k;
			g_rate_arr[i-1].frequency = g_rate_arr[i].frequency;
		}

		--g_rate_cnt;

		return;
	}

	transition_end_idx = g_rate_arr[rate_idx].transition_start_idx + g_rate_arr[rate_idx].transition_count; // last transition of same rate type

	// swap transition_end_idx into the position atom_idx:offset_idx occupied
	// ENHANCE: this is the same shit that happens when count==0
	g_transition_arr[transition_idx]->atom_idx = g_transition_arr[transition_end_idx]->atom_idx;
	g_transition_arr[transition_idx]->offset_idx = g_transition_arr[transition_end_idx]->offset_idx;

	if (transition_idx != transition_end_idx) 
		g_atom_arr[g_transition_arr[transition_idx]->atom_idx]->transition_indices[g_transition_arr[transition_idx]->offset_idx] = transition_idx;

	// shift all other transition lists

	for (i=rate_idx+1;i < g_rate_cnt;++i)
	{ // ENHANCE: again, looks like the same shit that happens when count==0
		--g_rate_arr[i].transition_start_idx;

		transition_idx = g_rate_arr[i].transition_start_idx;
		transition_end_idx = g_rate_arr[i].transition_start_idx+g_rate_arr[i].transition_count;

		g_transition_arr[transition_idx]->atom_idx = g_transition_arr[transition_end_idx]->atom_idx;
		g_transition_arr[transition_idx]->offset_idx = g_transition_arr[transition_end_idx]->offset_idx;

		g_atom_arr[g_transition_arr[transition_idx]->atom_idx]->transition_indices[g_transition_arr[transition_idx]->offset_idx] = transition_idx;
	}

	free(g_transition_arr[g_transition_cnt]);			// free up the very last member of the last g_transition_arr

	return;
}


/******************************************************************************/
/******************************************************************************/

void check_system(struct SimulationState *ss)
{
	int i, j, k, m, n, mm;
	int errors;

	int next_x, next_y, next_z;

	int nnx, nny, nnz;

	int xzone, yzone, zzone;

	double xx1[3], xx2[3], xx3[3];
	double yy1[3], yy2[3], yy3[3];

	// does a careful check to make sure that a system is ready to be simulated.
	// assumptions:  (1) all real atoms are in the places they think they are
	// (2) all buried atoms are actually buried.

	// first, remove all atoms from the transition list.  We'll add them after we check neighbors

	for (j=0;j<g_atom_cnt;++j)
		for (i=0;i<g_max_neighbors + 1;++i)				// extra 1 for evaporation
			{
				if (g_atom_arr[j]->transition_indices[i] != -1)		// something can happen in the "i" direction
					take_off_transition_list(j, i);
			}

	// for each atom, cycle through neighbor coordinates, and and reconcile occupancy

	for (j=0;j<g_atom_cnt;++j)
		for (i=0;i<g_max_neighbors;++i)
		{
   			next_x = g_atom_arr[j]->lattice[0] + jump_offset[i].dx;
			next_y = g_atom_arr[j]->lattice[1] + jump_offset[i].dy;
	        next_z = g_atom_arr[j]->lattice[2] + jump_offset[i].dz;

			adjust_pbc(&next_x, &next_y, &next_z);

			k = atom_at(next_x, next_y, next_z, ss);

			if (k >= 0)
			{
				// an atom has been found atom_idx this neighbor site.

				g_atom_arr[j]->occupied_neighbor_sites[i] = k;
				g_atom_arr[k]->occupied_neighbor_sites[opposite_offset[i]] = j;
			}
		}

	// now we'll reconcile buried atoms.
	// does this process need to happen if nothing is buried???
	do
	{
		errors = 0;

		for (j=0;j<g_atom_cnt;++j)
		for (i=0;i<g_max_neighbors;++i)
		{
			if (g_atom_arr[j]->occupied_neighbor_sites[i] == -2)
			{
				// find coordinate of buried atom

				next_x = g_atom_arr[j]->lattice[0] + jump_offset[i].dx;
				next_y = g_atom_arr[j]->lattice[1] + jump_offset[i].dy;
		        next_z = g_atom_arr[j]->lattice[2] + jump_offset[i].dz;

				adjust_pbc(&next_x, &next_y, &next_z);

				for (k=0;k<g_max_neighbors;++k)
				{
					nnx = next_x + jump_offset[k].dx;
					nny = next_y + jump_offset[k].dy;
					nnz = next_z + jump_offset[k].dz;

					adjust_pbc(&nnx, &nny, &nnz);

					m = atom_at(nnx, nny, nnz, ss);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = g_atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]];
						if (n != -2)
						{
							if (n == -1) g_atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]] = -2;
							if (n == -3) g_atom_arr[j]->occupied_neighbor_sites[i] = -3;		// random trumps buried

							++errors;
						}
					}
				}
			}
		}
	}
	while (errors != 0);

	// now we'll reconcile random buried atoms.

	do
	{
		errors = 0;

		for (j=0;j<g_atom_cnt;++j)
		for (i=0;i<g_max_neighbors;++i)
		{
			if (g_atom_arr[j]->occupied_neighbor_sites[i] == -3)
			{
				// find coordinate of buried atom

				next_x = g_atom_arr[j]->lattice[0] + jump_offset[i].dx;
				next_y = g_atom_arr[j]->lattice[1] + jump_offset[i].dy;
		        next_z = g_atom_arr[j]->lattice[2] + jump_offset[i].dz;

				adjust_pbc(&next_x, &next_y, &next_z);

				for (k=0;k<g_max_neighbors;++k)
				{
					nnx = next_x + jump_offset[k].dx;
					nny = next_y + jump_offset[k].dy;
					nnz = next_z + jump_offset[k].dz;

					adjust_pbc(&nnx, &nny, &nnz);

					m = atom_at(nnx, nny, nnz, ss);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = g_atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]];
						if (n != -3)
						{
							g_atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]] = -3;
							++errors;
						}
					
					}
				}
			}
		}
	}
	while (errors != 0);

	// now let's bury any atoms that should be buried - DON'T WANT THIS NOW!
	/*for (j=0;j<g_atom_cnt;++j)
	{
		k = 0;		// k will be the number of buried or occupied neighbors

		for (i=0;i<g_max_neighbors;++i)
		{
			if (atom[j]->occupied_neighbor_sites[i] >= 0)
			{	
				if (atom[j]->type == atom[atom[j]->occupied_neighbor_sites[i]]->type)
				++k;
			}

			if (atom[j]->occupied_neighbor_sites[i] == -2)
				++k;
		}

		// SPECIAL SC routine included here otherwise for second nearest neighbors?

		if (k == g_max_neighbors)
		{
			// bury atom j

			for (i=0;i<g_max_neighbors;++i)
			{
				m = atom[j]->occupied_neighbor_sites[i];
				if (m >= 0)
					atom[m]->occupied_neighbor_sites[opposite_offset[i]] = -2;
			}
				
			// remove the atom from the atom list

			m = atom[j]->next_atom;
			n = atom[j]->previous_atom;

			if (n == -1)
			{
				// this is first atom on this list, so make the zone point to
				// the next element in the list.  Note that if the zone had only
				// one element, i should be -1, which will alert the offset that
				// the zone is empty

				findzone (&xzone, &yzone, &zzone, atom[j]->lattice[0], atom[j]->lattice[1],atom[j]->lattice[2]);

				zone[xzone][yzone][zzone].offset = m;

				if (m != -1)
					atom[m]->previous_atom = -1;
			}
			else
			{
				if (m == -1)
					{
						// this is the last element on this list,
						atom[n]->next_atom = -1;
					}
				else
					{
						// atom is embedded in the list; nothing special needs be done
	
						atom[m]->previous_atom = n;
						atom[n]->next_atom = m;
					}
			}

			if (j != (g_atom_cnt-1))
					move_atom((g_atom_cnt-1), j, ss);

			free(atom[g_atom_cnt-1]);
			--g_atom_cnt;
			--j;
		}
	}*/

	// now we can recalculate diffusion rates

	for (j=0;j<g_atom_cnt;++j)
		refresh_transitions(j, ss);

	return;
}


/******************************************************************************/
/******************************************************************************/
// calculates surface diffusion rate constant
int calculate_surf_diffusion_rate(	int initial_configuration[],			// initial configuration array
									int final_configuration[],				// final configuration array
									int number_of_neighbors,				// number of 1st near neighbors in current xtal structure
									int atom_type,							// type of atom in consideration for transition
									double nnE[6],							// bond energy (type 2)-(type 2) [nearest neighbor energy]
									double temperature,						// system temperature
									double overpotential,					// system overpotential
									double *rate)							// return value - rate constant k
{ // ENHANCE: pass the number of nearest neighbors? (best practice)
	int i;
	double energy = 0.0;
	
	// [ ]: why are these doubles? presumably so they don't get implicitly cast as such when used in calculation
	// total number of neighbors in initial/final
	double neighbor_cnt_initial = 0.0;
	double neighbor_cnt_final = 0.0;
	
	// number of type A/B/C [index 012] neighboring atoms in initial/final configuration
	double neighbor_type_cnt_initial[] = {0.0, 0.0, 0.0};
	double neighbor_type_cnt_final[] = {0.0, 0.0, 0.0};
	
	int nnE_A_idxs[3] = {0, 1, 2}; // indices of A-A, A-B, A-C bonds in nnE array
	int nnE_B_idxs[3] = {1, 3, 4}; // indices of B-A, B-B, B-C bonds
	int nnE_C_idxs[3] = {2, 4, 5}; // indices of C-A, C-B, C-C bonds

	// [ ]: why is this static, not const? it should also be a simulation input
	static double b_anisotropy_factor[12] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; // "optional" anistropy factor, indices match jump_offset
	
	int neighbor_type; //type of atom for nearest neighbor
	// XXX: commented prints
	//printf("before teh switch: atom type = %d\n", atom_type);

	/*printf("final configuration: ");
	for (i=0;i<g_max_neighbors;++i)
		printf("%d ", final_configuration[i]);
	printf("\n");*/

	// ENHANCE: lol these are identical except for nnE[nne*_index...], do it better (use fxn)
	switch(atom_type)
	{
		case 1:
			for (i=0; i<number_of_neighbors; ++i)
			{
				neighbor_type = initial_configuration[i];
				//printf("init neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++neighbor_type_cnt_initial[neighbor_type - 1]; //number of type A/B/C neighboring atoms in init configuration
					energy += nnE[nnE_A_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //A-A/B/C bond - grab correct index of nnE array from nne*_index array based on neighbor type
				}
				neighbor_type = final_configuration[i];
				//printf("final neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++neighbor_type_cnt_final[neighbor_type - 1]; //number of type A/B/C neighboring atoms in final configuration
					//energy -= nnE[nnE_A_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //A-A/B/C bond, told to leave in and comment out
				}
			}
			break;


		case 2:

			for (i=0;i<number_of_neighbors;++i)
			{
				neighbor_type = initial_configuration[i];
				//printf("init neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++neighbor_type_cnt_initial[neighbor_type - 1]; //number of type A/B/C neighboring atoms in init configuration
					energy += nnE[nnE_B_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //B-A/B/C bond
				}
				neighbor_type = final_configuration[i];
				//printf("final neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++neighbor_type_cnt_final[neighbor_type - 1]; //number of type A/B/C neighboring atoms in final configuration
					//energy -= nnE[nnE_B_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //B-A/B/C bond, told to leave in and comment out
				}
			}
			break;

		case 3:
			for (i=0;i<number_of_neighbors;++i)
			{
				neighbor_type = initial_configuration[i];
				//printf("init neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++neighbor_type_cnt_initial[neighbor_type - 1]; //number of type A/B/C neighboring atoms in init configuration
					energy += nnE[nnE_C_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //C-A/B/C bond
				}
				neighbor_type = final_configuration[i];
				//printf("final neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++neighbor_type_cnt_final[neighbor_type - 1]; //number of type A/B/C neighboring atoms in final configuration
					//energy -= nnE[nnE_C_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //C-A/B/C bond, told to leave in and comment out
				}
			}
			break;

	}

	//printf("after the switch\n");
	for (i=0;i<3;++i)
	{
		neighbor_cnt_initial += neighbor_type_cnt_initial[i];
		neighbor_cnt_final += neighbor_type_cnt_final[i];
	}
	// these override the previous energy sum
	if (neighbor_cnt_initial == 0)
	{
		// no neighbors - this condition corresponds to a diffuser walking through a lattice (a lattice gas)
		
		energy = -1.0;

	}

	if ((neighbor_cnt_initial > 0) && (neighbor_cnt_final <= 1.0))
	{ // [ ]: this is preventing evaporation (when no final neighbors)
		energy = 1000.;				// final configuration has no near neighbors,
									// so this effectively corresponds to an evaporation-like event.
									// Don't let it happen!
	}

	// ENHANCE: replace calculating the exp with memoizing up the value (uhash?) -> speedup?
	// BUG: why was overpotential removed? Isn't it necessary for linear sweep stuff?
	//*rate = 1e13*exp(-energy/(kBoltz*ss->temperature)) //+ 1e-4*exp(-(energy-ss->overpotential)/(kBoltz*ss->temperature));
	*rate = 1e13*exp(-energy/(kBoltz*temperature));
	//printf("rate = %le\n", *rate);
	return 0;
}

/******************************************************************************/
/******************************************************************************/

int calculate_evaporation_rate(	int initial_configuration[],			// initial configuration array of atom's nearest neighbors
									int number_of_neighbors,				// number of 1st near neighbors in current xtal structure
									int atom_type,							// type of atom in consideration for transition
									double nnE[6],							// bond energy (type 2)-(type 2)
									double temperature,						// system temperature
									double overpotential,					// system overpotential
									double *rate)							// return value
{
	int i;
	double energy;
	// XXX: not used
	double EAu = .5;
	double nscale = 0.0;

	static double b_anisotropy_factor[12] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};			// optional anistropy factor

	energy = 0.0;

	/*if ((ncsk == 1)||(g_evaporation_flag == false)) //I don't think this happens
	{
		energy = 1000.;
		*rate = 1e4*exp(-(energy)/(kBoltz*ss->temperature));
		return 1;
	}*/

	int nnE_A_idxs[3] = {0, 1, 2}; //indices of A-A, A-B, A-C bonds in nnE array
	int nnE_B_idxs[3] = {1, 3, 4}; //indices of B-A, B-B, B-C bonds
	int nnE_C_idxs[3] = {2, 4, 5}; //indices of C-A, C-B, C-C bonds
	int neighbor_type; //type of atom for nearest neighbor

	// ENHANCE: lol these are identical except for nnE[nne*_index...], be better (fxn)
	switch(atom_type)
	{ // [ ]: how much of this is duplicated with calculate_surf_diffusion_rate
		case 1:
			if (g_solubility[0]) {
				//A can evaporate
				for (i=0;i<number_of_neighbors;++i)
				{ // calculate energy of initial state before evaporation
					neighbor_type = initial_configuration[i];
					if (neighbor_type > 0)
						energy += nnE[nnE_A_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //bonding with A
				}
			}
			else
				energy = 1000.; //A cannot evaporate
			*rate = 1e4*exp(-(energy-overpotential)/(kBoltz*temperature)); // [ ]: dissolution/evaporation equation
			break;

		case 2:
			if (g_solubility[1]) {
				//B can evaporate
				for (i=0;i<number_of_neighbors;++i)
				{
					neighbor_type = initial_configuration[i];
					if (neighbor_type > 0)
						energy += nnE[nnE_B_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //bonding with B
				}
			}
			else
				energy = 1000.; //B cannot evaporate
			*rate = 1e4*exp(-(energy-overpotential)/(kBoltz*temperature));
			break;

		case 3:
			if (g_solubility[2]) {
				//C can evaporate
				for (i=0;i<number_of_neighbors;++i)
				{
					neighbor_type = initial_configuration[i];
					if (neighbor_type > 0)
						energy += nnE[nnE_C_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //bonding with C
				}
			}
			else
				energy = 1000.; //C cannot evaporate
			*rate = 1e4*exp(-(energy-overpotential)/(kBoltz*temperature));
			break;
	}
		
	return 0;
}
