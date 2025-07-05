#include "Simulation.h"
#include "Random.h"
#include "Simulation_Aux.h"
#include "FileIO.h"
#include "Atoms.h"
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

int rate_skip; // transition list binary search index
int adatom_before; // XXX: never used?

int lastxt, lastyt, lastzt; // containers for coordinates of a next step
bool checkpoint_reached = false;

// ENHANCE: pass struct with all simulation parameters as argument
unsigned long perform_simulation(struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls) //potentially FILE* as arguments
{
	//printf("Simulation starting\n");
	long int i;
	int j, k;

	double cur_stime, prev_stime; // prev and current times, for overpotential moving

	double transition_type_probability;

	//double vap;

	int which_one;
	int atom_number, jump_vector;
	int natn; //changed from nan bc keyword
	int atype;

	int moved_flag = true;

	int framenum = 0; // counter/id for number of outputs / output files

	ss->elapsed_stime = 0.0;
	//writes data time intervals and run time in original code

	//simulation_is_going = true;

	//print system and zone sizes to the file?

	/*if (num_sims > 0) //this probably will not happen
	{
		do_initialize_simulation(simulation_type);
		elapsed_stime = 0.0;
	}*/ 

	// initialize simulation kinetics, and draw a picture

	//printf("transition time!\n");
	for (i=0; i < ss->atom_cnt; ++i)	
		// resets all kinetic paramters
		refresh_transitions(i, ss, se);
	// TODO: move this out from here, to only where output is necessary
	organize(ss->atom_arr, ss->atom_cnt); //replacement for copy_xyz_to_coord but might not be necessary
	
	ss->total_atoms_dissolved = 0;

	// initial state
	calculate_internal_energy(ss->atom_arr, ss->atom_cnt, &ss->total_internal_energy, se);
	
	output_log_file(ls->sim_log_file, framenum, ss->elapsed_stime, ss->temperature, ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
	write_xyz_file(ss, ls->position_log_prefix, framenum);
	framenum++;

	unsigned long int iter = 1; // iteration count

	prev_stime = 0.0; //needs to happen outside of loop
	//printf("about to start the loop\n");
	
	bool simulation_end = false;

	while (!simulation_end)
	{
		if (iter % 100 == 0)
			printf("iteration %ld, time %lf\n", iter, ss->elapsed_stime);
		if (ss->simulation_should_kill_itself) // abort simulation (only happens if atoms overlap)
		{
			//find_average_curvature(); // XXX: no longer valid

			calculate_internal_energy(ss->atom_arr, ss->atom_cnt, &ss->total_internal_energy, se);
			output_log_file(ls->sim_log_file, framenum, ss->elapsed_stime, ss->temperature, ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
			write_xyz_file(ss, ls->position_log_prefix, framenum);

			ss->simulation_should_kill_itself = false;

			// TODO: remove from here
			organize(ss->atom_arr, ss->atom_cnt); // XXX: remove from here

			return 1; //return 1 b/c error
		}
		
		// TODO: only do this for the atom that was changed and its neighbors
		for (j=0; j < ss->atom_cnt; ++j)
			refresh_transitions(j, ss, se); // resets all kinetic paramters
		//does this need to happen?
		
		// increment the elapsed time
		compute_transition_array(ss);
		ss->elapsed_stime -= log(drandj(&rand_seed)) / ss->frequency_sum;

		// if (ss->elapsed_stime >= ss->run_stime) // simulation has gone past time
			// break; //get outta here before I make a new transition

		// pick the type of transition to occur
		transition_type_probability = drandj(&rand_seed);
		rate_skip = ss->rate_cnt / 2;
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
				which_one = ss->rate_arr[k].transition_start_idx + (int)(drandj(&rand_seed) * (double)ss->rate_arr[k].transition_count);

				// which_one gives the location of the transition_arr, which gives
				// the info about the specific atom
				atom_number = ss->transition_arr[which_one]->atom_idx;
				jump_vector = ss->transition_arr[which_one]->offset_idx;

				// if jump_vector == se->max_neighbors then the atom is going to evaporate
				moved_flag = true;

				adatom_before = 0;

				if (jump_vector != se->max_neighbors) //diffusion
				{
					//printf("the transition is diffusion of atom %d, jumping from %lf %lf %lf ", atom_number, atom[atom_number]->lattice[0], atom[atom_number]->lattice[1], atom[atom_number]->lattice[2]);
					// coordinates atom is jumping to
					lastxt = ss->atom_arr[atom_number]->lattice[0] + jump_offset[jump_vector].dx;
					lastyt = ss->atom_arr[atom_number]->lattice[1] + jump_offset[jump_vector].dy;
					lastzt = ss->atom_arr[atom_number]->lattice[2] + jump_offset[jump_vector].dz;

					adjust_pbc(&lastxt, &lastyt, &lastzt, se);

					atype = ss->atom_arr[atom_number]->type;

					// moves atom?
					remove_atom(atom_number, ss, se);
					natn = add_atom(lastxt, lastyt, lastzt, atype, NORMAL, ss, se);
					//printf("and jumping to %lf %lf %lf\n", lastxt, lastyt, lastzt);
				}
				else				// dissolution
				{
					if (se->solubility[ss->atom_arr[atom_number]->type - 1] == true)	// atoms are dissolved based on input specs!
					{
						++ss->total_atoms_dissolved;
						remove_atom(atom_number, ss, se);	// evaporate the atom
						//printf("the transition was dissolution of atom %d\n", atom_number);
					}
				}
			}
			else if (transition_type_probability < (ss->transition_probability.lbound[j]))
			{
				//search to the left
				rate_skip = rate_skip/2;
				if (rate_skip == 0) rate_skip = 1;
				j -= rate_skip;
			}
			else if (transition_type_probability >= (ss->transition_probability.ubound[j]))
			{
				//search to the right
				rate_skip = rate_skip/2;
				if (rate_skip == 0)
					rate_skip = 1;
				j += rate_skip;
				if (j == ss->rate_cnt) //no more options!
					break;
			}
		}

		if (moved_flag == false) //only happens iff jump_vector == se->max_neighbors
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

					rrp[0] = ssr*cos(ta2)*cos(ta1);
					rrp[1] = ssr*cos(ta2)*sin(ta1);
					rrp[2] = ssr*sin(ta2);

					// now invert xr, yr, zr into lattice vectors;

					vecmul(rrp, invert_primitive_basis, rw);

					rwx = (int)rw[0] + ssx/2;
					rwy = (int)rw[1] + ssy/2;
					rwz = (int)rw[2] + ssz/2;

				} while (atom_at(rwx, rwy, rwz) >= 0);

				//if (atom_at(rwx, rwy, rwz) >= 0) goto newt; (made redundant with do while)

				natn = add_atom(rwx, rwy, rwz, atype, NORMAL);
			}*/
		}
		 
		// after iteration, log if necessary
		// TODO: implement the checkpoint lists
		if ((ls->analysis_type == REGULAR_TIME_INTERVALS) || (ls->analysis_type == LN_TIME_INTERVALS))
			checkpoint_reached = (ss->elapsed_stime >= ls->next_log_checkpoint);
		else if (ls->analysis_type == ITERATION_INTERVALS)
			checkpoint_reached = (iter >= ls->next_log_checkpoint);

		if (checkpoint_reached)
		{
			organize(ss->atom_arr, ss->atom_cnt); //replaced but do i really need it
			
			//record the elapsed time in a file here
			printf("writing file %d: elapsed_stime = %lf\n", framenum, ss->elapsed_stime);
			
			calculate_internal_energy(ss->atom_arr, ss->atom_cnt, &ss->total_internal_energy, se);
			output_log_file(ls->sim_log_file, framenum, ss->elapsed_stime, ss->temperature, ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
			write_xyz_file(ss, ls->position_log_prefix, framenum);
			
			if (ls->analysis_type == REGULAR_TIME_INTERVALS)
			{
				// bring next_log_checkpoint up to and one step beyond elapsed_stime
				while (ls->next_log_checkpoint <= ss->elapsed_stime)
					ls->next_log_checkpoint += ls->log_interval;

			}
			else if (ls->analysis_type == LN_TIME_INTERVALS)
			{
				while (ls->next_log_checkpoint <= ss->elapsed_stime)
					ls->next_log_checkpoint *= ls->log_interval;
			}
			else if (ls->analysis_type == ITERATION_INTERVALS)
			{
				ls->next_log_checkpoint += ls->log_interval;
			}

			if (se->overpotential_ramp_rate != 0.0)
			{
				cur_stime = ss->elapsed_stime;
				ss->overpotential += (cur_stime-prev_stime) * se->overpotential_ramp_rate;
				prev_stime = ss->elapsed_stime;
				for (j=0; j < ss->atom_cnt; ++j)	
					refresh_transitions(j, ss, se); // resets all kinetic paramters

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
		fprintf(ls->sim_log_file, "reached final iteration and terminated\n");
	
	//write elapsed_stime to mark finish
	//TODO: finish IO
	calculate_internal_energy(ss->atom_arr, ss->atom_cnt, &ss->total_internal_energy, se);
	output_log_file(ls->sim_log_file, framenum, ss->elapsed_stime, ss->temperature, ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
	write_xyz_file(ss, ls->position_log_prefix, framenum);

	printf("Finished simulation\n"); //move this to the log file
		
	//simulation_is_going = false;
	return 0;
}

/******************************************************************************/
/******************************************************************************/
// updates ss->transition_probability (weighted rate list, used to choose event)
void compute_transition_array(struct SimulationState* ss)
{
	int i;
	int x;

	int total_lists = 0;

	double sum_of_rate_populations = 0.0;
	ss->frequency_sum = 0.0;

	for (i=0; i < ss->rate_cnt; ++i)
	{
		if (ss->rate_arr[i].transition_count != 0)
		{
			ss->rate_arr[i].frequency = ss->rate_arr[i].k*(double)ss->rate_arr[i].transition_count;
			ss->frequency_sum += ss->rate_arr[i].frequency;
			sum_of_rate_populations += ss->rate_arr[i].transition_count;

			ss->transition_probability.listnum[total_lists] = i;
			++total_lists;
		}
	}

	/*if (deposition_type == DEPOSITION_TYPE_RAINFALL)
		frequency_sum += (double)ssx*(double)ssy*(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c);
	else if (deposition_type == DEPOSITION_TYPE_RANDOM_WALKER)
		frequency_sum += (4*PI*ssr*ssr)*(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c);*/ //deposition is vestigial

	// now compute bounds for jump probabilities

	double current_probability = 0.0;

	for (i=0;i<total_lists;++i)
	{
		ss->transition_probability.lbound[i] = current_probability;
		x = ss->transition_probability.listnum[i];
		current_probability += ss->rate_arr[x].frequency / ss->frequency_sum;
		ss->transition_probability.ubound[i] = current_probability;
	}

	return;
}

int get_bond_index(int a, int b, struct SimulationEnv* se)
{
	// aa, ab, ac; bb, bc; cc [num_elements=3]
	// 00, 01, 02; 11, 12; 22
	// assume 0-indexed
	int first, second;

	if (a > b)
	{
		first = a;
		second = b;
	}
	else {
		first = b;
		second = a;
	}
	
	// a=1, b=2 -> (1*3)+(2-1)=4
	return (a * (se->num_elements)) + (b-a);
}

/******************************************************************************/
/******************************************************************************/
// updates [atom_arr[atom_idx], transition_arr[i], rate_arr[i].transition_start_idx], initializes rate_arr
int refresh_transitions(int atom_idx, struct SimulationState* ss, struct SimulationEnv* se) // atom_idx = index on atom list
{
	int i, j;
	int rate_idx; // position of rate rate in rate list rate_arr[]
	double rate, evap_rate; // rate constant from bond-breaking model
	int next_x, next_y, next_z;
	Atom *cur_atom = ss->atom_arr[atom_idx];

	int atom_rates_cnt;	// this is returned as the number of transitions (jump_offsets) this atom can undergo, excluding evaporation
	int start_config[MAXIMUM_NUMBER_OF_NEIGHBORS]; // -1 if empty, type if filled
	int end_config[MAXIMUM_NUMBER_OF_NEIGHBORS];
	
	unsigned char *atom_env = (unsigned char*) calloc(se->num_nn_levels * se->num_bond_types, sizeof(unsigned char));
	// unsigned char env_hash[se->num_nn_levels * se->num_bond_types];

	// first, remove all mention of this atom from transition list
	for (i=0; i < se->max_neighbors + se->dissolution; ++i) // extra 1 for evaporation
	{
		if (ss->atom_arr[atom_idx]->transition_indices[i] != -1) // something can happen in the "i" direction
			take_off_transition_list(atom_idx, i, ss);
	}

	// cycle through neighbor coordinates, and check if there is an atom there
	for (i=0;i < se->max_neighbors;++i)
	{
		ss->atom_arr[atom_idx]->occupied_neighbor_sites[i] = 0;
   		next_x = ss->atom_arr[atom_idx]->lattice[0] + jump_offset[i].dx;
		next_y = ss->atom_arr[atom_idx]->lattice[1] + jump_offset[i].dy;
        next_z = ss->atom_arr[atom_idx]->lattice[2] + jump_offset[i].dz;

		adjust_pbc(&next_x, &next_y, &next_z, se);

        j = atom_at(next_x, next_y, next_z, ss->atom_arr, ss->zone_arr, se);

		// update atom_env
		if (j >= 0)
		{
			ss->atom_arr[atom_idx]->occupied_neighbor_sites[i] = j;
			
			int bond_idx = get_bond_index(ss->atom_arr[atom_idx]->type, ss->atom_arr[j]->type, se);
			int env_idx = get_env_index(1, bond_idx, se);
			// 0 for 1st nn
			atom_env[env_idx]++;
		}

		// // if no atom atom_idx that position but occ_neighbor array says there is, fix it
		// if ((j == -1) && (ss->atom_arr[atom_idx]->occupied_neighbor_sites[i] >= 0))
		// 	ss->atom_arr[atom_idx]->occupied_neighbor_sites[i] = -1;
	}

	// cycle through the neighbor sites.  if there's an empty one, calculate the transition rate to it

	atom_rates_cnt = 0;

	int intial_config_neighbor_cnt = get_initial_configuration2(atom_idx, 0, se->max_neighbors, ss->atom_arr, start_config);		// k is number of near neighbors

	if (intial_config_neighbor_cnt == se->max_neighbors) // skip calculating a rate of a fully coordinated atom
		return atom_rates_cnt;
	
	for (i=0; i < se->max_neighbors; ++i) // create transitions to each unoccupied neighbor
	{
		if (ss->atom_arr[atom_idx]->occupied_neighbor_sites[i] == -1) // if unoccupied, consider transition
		{ 
			// end_config is only used to identify the transitions that are functionally evaporations (no neighbors in end configuration)
			int final_config_neighbor_cnt = get_final_configuration2(atom_idx, i, ss, se, end_config);
			//printf("surf diffusion\n");
			calculate_surf_diffusion_rate(		
				start_config, // ENHANCE: make this look prettier
				end_config,
				ss->atom_arr[atom_idx]->type,
				ss->temperature,
				ss->overpotential,
				&rate,
				se
			);
					
			++atom_rates_cnt;
			//printf("adding step\n");
			if ((rate_idx = is_on_rate_list(ss, rate)) != -1)
			{ // if rate constant already in rate list, use that rate list index
				add_to_transition_list(rate_idx, atom_idx, i, ss, se);
			}
			else
			{
				// the transition rate with rate constant rate turned out to be a new one
				// add to rate list/array
				rate_idx = create_new_transition(ss, rate);
				add_to_transition_list(rate_idx, atom_idx, i, ss, se);
			}
		}
	}	

	//printf("calcualte evap\n");
	calculate_evaporation_rate(	
		start_config,
		ss->atom_arr[atom_idx]->type,
		ss->temperature,
		ss->overpotential,
		&rate,
		se
	);

	//replace with the lines below because it's more obvious
	/*if ((rate_idx = is_on_rate_list(rate)) != -1)
	{
		add_to_transition_list(rate_idx, atom_idx, se->max_neighbors);
	}
	else  // the transition rate to that spot turned out to be a new one.
	{		
		rate_idx = create_new_transition(rate);
		add_to_transition_list(rate_idx, atom_idx, se->max_neighbors);
	}*/

	//printf("end stuff\n");
	rate_idx = is_on_rate_list(ss, rate);

	if (rate_idx == -1)
		rate_idx = create_new_transition(ss, rate); //the transition rate to the spot is a new one!

	add_to_transition_list(rate_idx, atom_idx, se->max_neighbors, ss, se); // evaporation is considered to be last in jump_offset (not really in array but uses that index number)
	//printf("gonna return %d\n", atom_rates_cnt);
	return atom_rates_cnt;						// gives number of current transitions for that atom
}

/******************************************************************************/
/******************************************************************************/

//is there a better way of doing this?
// checks if the rate constant is already in the rate list
int is_on_rate_list(struct SimulationState* ss, double rate) 
{
	int i;

	for (i=0; i < ss->rate_cnt; ++i)
		if (rate == ss->rate_arr[i].k) return i;

	return -1;
}

/******************************************************************************/
/******************************************************************************/
// create new Rate struct in rate array
// updates rate_array[rate_cnt], rate_cnt
int create_new_transition(struct SimulationState* ss, double rate)
{
	ss->rate_arr[ss->rate_cnt].k = rate;
	ss->rate_arr[ss->rate_cnt].transition_start_idx = ss->transition_cnt;
	ss->rate_arr[ss->rate_cnt].transition_count = 0;

	++ss->rate_cnt;

	return (ss->rate_cnt-1);
}

/******************************************************************************/
/******************************************************************************/

// add to rate_arr[rate_idx] the atom atom_idx going in direction offset_idx
// updates transition_arr, rate_arr[rate_idx].transition_count, atom_arr[atom_idx]->transition_indices[offset_idx]
void add_to_transition_list(int rate_idx, int atom_idx, int offset_idx, struct SimulationState* ss, struct SimulationEnv* se) // rate_arr index, atom_arr index, jump_offset index
{
	int i; // loop variable
	int n;
	int initial_transition_index, final_transition_index; // initial and final transition_arr index

	// make room for the new arrival

	ss->transition_arr[ss->transition_cnt] = (Transition *)malloc(sizeof(Transition));	// adds entry to the end of the list
	if (ss->transition_arr[ss->transition_cnt] == NULL)
	{
		// TODO: free mallocs before exiting
		fprintf(stderr, "Couldn't allocate memory for atom %d: %s\n", ss->transition_cnt, strerror(errno));
		exit(errno);
	}
	if ((unsigned int) ss->transition_cnt > se->max_transitions)
	{
		fprintf(stderr, "More transitions (%d) than allocated in transition array (%llu)\n", ss->transition_cnt, se->max_transitions);
		exit(errno);
	}

	// what is this
	//		final_transition_index = rate_arr[rate_cnt-1].transition_start_idx + rate_arr[rate_cnt-1].number;	
	//		transition_arr[final_transition_index] = (Transition *)malloc(sizeof(Transition));

	for (i = ss->rate_cnt-1; i > rate_idx; --i)
	{ // [ ]: what does this do? is this the same as in remove_transition?
		initial_transition_index = ss->rate_arr[i].transition_start_idx;
		final_transition_index = initial_transition_index + ss->rate_arr[i].transition_count;

		ss->transition_arr[final_transition_index]->atom_idx = ss->transition_arr[initial_transition_index]->atom_idx;
		ss->transition_arr[final_transition_index]->offset_idx = ss->transition_arr[initial_transition_index]->offset_idx;

		if (initial_transition_index != final_transition_index)
			ss->atom_arr[ss->transition_arr[final_transition_index]->atom_idx]->transition_indices[ss->transition_arr[final_transition_index]->offset_idx] = final_transition_index;

		++ss->rate_arr[i].transition_start_idx;
	}

	// add new arrival

	n = ss->rate_arr[rate_idx].transition_start_idx + ss->rate_arr[rate_idx].transition_count;

	++ss->rate_arr[rate_idx].transition_count;

	ss->transition_arr[n]->atom_idx = atom_idx;
	ss->transition_arr[n]->offset_idx = offset_idx;

	ss->atom_arr[atom_idx]->transition_indices[offset_idx] = n;

	++ss->transition_cnt;

	return;
}

/******************************************************************************/
/******************************************************************************/
// updates atom_arr[atom_idx], transition_arr[i], rate_arr[i].transition_start_idx
void take_off_transition_list(int atom_idx, int offset_idx, struct SimulationState* ss)	// removes atom jumping in the jump_offset[offset_idx] direction
{
	int i;
	int rate_idx, transition_idx, transition_end_idx;

	// find out what Rate in rate_arr this is

	transition_idx = ss->atom_arr[atom_idx]->transition_indices[offset_idx]; // old position on transition list, to be removed

	for (i=0; i < ss->rate_cnt; ++i)
		if (transition_idx < (ss->rate_arr[i].transition_start_idx + ss->rate_arr[i].transition_count))
		{
			rate_idx = i;
			break;
		}

	// remind atom it can no longer jump

	ss->atom_arr[atom_idx]->transition_indices[offset_idx] = -1;
	// TODO: do these later, after it's been removed and transition_arr has been rearranged
	--ss->transition_cnt;

	// rate_idx points to the current rate list it's on.  decrement the number of atoms in that list
	// and clean up.  If the list is empty, remove it.

	--ss->rate_arr[rate_idx].transition_count;
	// [ ]: wtf is going on here 
	if (ss->rate_arr[rate_idx].transition_count == 0) // if list is empty
	{
		for (i = rate_idx + 1; i < ss->rate_cnt; ++i)
		{
			--ss->rate_arr[i].transition_start_idx; // move rate_arr offsets of larger indicies down one
			// make transition at new start index (which is of different Rate than old start index) have same atom and offset as new end index (which is of same Rate as old end index)
			transition_idx = ss->rate_arr[i].transition_start_idx;
			transition_end_idx = ss->rate_arr[i].transition_start_idx + ss->rate_arr[i].transition_count;		// count is always at least 1

			ss->transition_arr[transition_idx]->atom_idx = ss->transition_arr[transition_end_idx]->atom_idx;
			ss->transition_arr[transition_idx]->offset_idx = ss->transition_arr[transition_end_idx]->offset_idx;
			// update the transition index in the corresponding atom in atom_arr to have the new (lower) transition index 
			ss->atom_arr[ss->transition_arr[transition_idx]->atom_idx]->transition_indices[ss->transition_arr[transition_idx]->offset_idx] = transition_idx;
		}

		free(ss->transition_arr[ss->transition_cnt]);			// free up the very last member of the last transition_arr

		for (i = rate_idx + 1; i < ss->rate_cnt; ++i)
		{
			ss->rate_arr[i-1].transition_start_idx = ss->rate_arr[i].transition_start_idx;
			ss->rate_arr[i-1].transition_count = ss->rate_arr[i].transition_count;
			ss->rate_arr[i-1].k = ss->rate_arr[i].k;
			ss->rate_arr[i-1].frequency = ss->rate_arr[i].frequency;
		}

		--ss->rate_cnt;

		return;
	}

	transition_end_idx = ss->rate_arr[rate_idx].transition_start_idx + ss->rate_arr[rate_idx].transition_count; // last transition of same rate type

	// swap transition_end_idx into the position atom_idx:offset_idx occupied
	// ENHANCE: this is the same shit that happens when count==0
	ss->transition_arr[transition_idx]->atom_idx = ss->transition_arr[transition_end_idx]->atom_idx;
	ss->transition_arr[transition_idx]->offset_idx = ss->transition_arr[transition_end_idx]->offset_idx;

	if (transition_idx != transition_end_idx) 
		ss->atom_arr[ss->transition_arr[transition_idx]->atom_idx]->transition_indices[ss->transition_arr[transition_idx]->offset_idx] = transition_idx;

	// shift all other transition lists

	for (i = rate_idx + 1; i < ss->rate_cnt;++i)
	{ // ENHANCE: again, looks like the same shit that happens when count==0
		--ss->rate_arr[i].transition_start_idx;

		transition_idx = ss->rate_arr[i].transition_start_idx;
		transition_end_idx = ss->rate_arr[i].transition_start_idx + ss->rate_arr[i].transition_count;

		ss->transition_arr[transition_idx]->atom_idx = ss->transition_arr[transition_end_idx]->atom_idx;
		ss->transition_arr[transition_idx]->offset_idx = ss->transition_arr[transition_end_idx]->offset_idx;

		ss->atom_arr[ss->transition_arr[transition_idx]->atom_idx]->transition_indices[ss->transition_arr[transition_idx]->offset_idx] = transition_idx;
	}

	free(ss->transition_arr[ss->transition_cnt]);			// free up the very last member of the last transition_arr

	return;
}


/******************************************************************************/
/******************************************************************************/

void check_system(struct SimulationState* ss, struct SimulationEnv* se)
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

	for (j = 0; j < ss->atom_cnt; ++j)
		for (i=0;i<se->max_neighbors + 1;++i)				// extra 1 for evaporation
			{
				if (ss->atom_arr[j]->transition_indices[i] != -1)		// something can happen in the "i" direction
					take_off_transition_list(j, i, ss);
			}

	// for each atom, cycle through neighbor coordinates, and and reconcile occupancy

	for (j = 0; j < ss->atom_cnt; ++j)
		for (i=0; i < se->max_neighbors;++i)
		{
   			next_x = ss->atom_arr[j]->lattice[0] + jump_offset[i].dx;
			next_y = ss->atom_arr[j]->lattice[1] + jump_offset[i].dy;
	        next_z = ss->atom_arr[j]->lattice[2] + jump_offset[i].dz;

			adjust_pbc(&next_x, &next_y, &next_z, se);

			k = atom_at(next_x, next_y, next_z, ss->atom_arr, ss->zone_arr, se);

			if (k >= 0)
			{
				// an atom has been found atom_idx this neighbor site.

				ss->atom_arr[j]->occupied_neighbor_sites[i] = k;
				ss->atom_arr[k]->occupied_neighbor_sites[opposite_offset[i]] = j;
			}
		}

	// now we'll reconcile buried atoms.
	// does this process need to happen if nothing is buried???
	do
	{
		errors = 0;

		for (j = 0; j < ss->atom_cnt; ++j)
		for (i=0;i<se->max_neighbors;++i)
		{
			if (ss->atom_arr[j]->occupied_neighbor_sites[i] == -2)
			{
				// find coordinate of buried atom

				next_x = ss->atom_arr[j]->lattice[0] + jump_offset[i].dx;
				next_y = ss->atom_arr[j]->lattice[1] + jump_offset[i].dy;
		        next_z = ss->atom_arr[j]->lattice[2] + jump_offset[i].dz;

				adjust_pbc(&next_x, &next_y, &next_z, se);

				for (k=0;k<se->max_neighbors;++k)
				{
					nnx = next_x + jump_offset[k].dx;
					nny = next_y + jump_offset[k].dy;
					nnz = next_z + jump_offset[k].dz;

					adjust_pbc(&nnx, &nny, &nnz, se);

					m = atom_at(nnx, nny, nnz, ss->atom_arr, ss->zone_arr, se);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = ss->atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]];
						if (n != -2)
						{
							if (n == -1) 
								ss->atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]] = -2;
							if (n == -3)
								ss->atom_arr[j]->occupied_neighbor_sites[i] = -3;		// random trumps buried

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

		for (j = 0; j < ss->atom_cnt; ++j)
		for (i = 0; i < se->max_neighbors; ++i)
		{
			if (ss->atom_arr[j]->occupied_neighbor_sites[i] == -3)
			{
				// find coordinate of buried atom

				next_x = ss->atom_arr[j]->lattice[0] + jump_offset[i].dx;
				next_y = ss->atom_arr[j]->lattice[1] + jump_offset[i].dy;
		        next_z = ss->atom_arr[j]->lattice[2] + jump_offset[i].dz;

				adjust_pbc(&next_x, &next_y, &next_z, se);

				for (k=0;k<se->max_neighbors;++k)
				{
					nnx = next_x + jump_offset[k].dx;
					nny = next_y + jump_offset[k].dy;
					nnz = next_z + jump_offset[k].dz;

					adjust_pbc(&nnx, &nny, &nnz, se);

					m = atom_at(nnx, nny, nnz, ss->atom_arr, ss->zone_arr, se);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = ss->atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]];
						if (n != -3)
						{
							ss->atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]] = -3;
							++errors;
						}
					
					}
				}
			}
		}
	}
	while (errors != 0);

	// now let's bury any atoms that should be buried - DON'T WANT THIS NOW!
	/*for (j=0;j<atom_cnt;++j)
	{
		k = 0;		// k will be the number of buried or occupied neighbors

		for (i=0;i<se->max_neighbors;++i)
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

		if (k == se->max_neighbors)
		{
			// bury atom j

			for (i=0;i<se->max_neighbors;++i)
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

			if (j != (atom_cnt-1))
					move_atom((atom_cnt-1), j);

			free(atom[atom_cnt-1]);
			--atom_cnt;
			--j;
		}
	}*/

	// now we can recalculate diffusion rates

	for (j = 0; j < ss->atom_cnt; ++j)
		refresh_transitions(j, ss, se);

	return;
}


/******************************************************************************/
/******************************************************************************/
// calculates surface diffusion rate constant
int calculate_surf_diffusion_rate(	
	int initial_configuration[],			// initial configuration array
	int final_configuration[],				// final configuration array
	int atom_type,							// type of atom in consideration for transition
	double temperature,						// system temperature
	double overpotential,					// system overpotential
	double *rate,							// return value - rate constant k
	struct SimulationEnv* se
)
{ // ENHANCE: pass the number of nearest neighbors? (best practice)
	int i;
	double energy = 0.0;
	
	// [ ]: why are these doubles? presumably so they don't get implicitly cast as such when used in calculation
	// total number of neighbors in initial/final
	int neighbor_cnt_initial = 0;
	int neighbor_cnt_final = 0;
	
	// number of type A/B/C [index 012] neighboring atoms in initial/final configuration
	int neighbor_type_cnt_initial[] = {0, 0, 0};
	int neighbor_type_cnt_final[] = {0, 0, 0};

	// [ ]: why is this static, not const? it should also be a simulation input
	static double b_anisotropy_factor[12] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; // "optional" anistropy factor, indices match jump_offset
	
	int neighbor_type; //type of atom for nearest neighbor

	int bond_idx, env_idx;
	for (i=0; i<se->max_neighbors; ++i)
	{
		bond_idx = get_bond_index(atom_type, neighbor_type, se);
		env_idx = get_env_index(1, bond_idx, se); // TODO: hard-coded 1st nn

		neighbor_type = initial_configuration[i];
		if (neighbor_type > 0) {
			++neighbor_type_cnt_initial[neighbor_type - 1]; //number of type A/B/C neighboring atoms in init configuration
			energy += se->nnEa[env_idx] * b_anisotropy_factor[i]; //A-A/B/C bond - grab correct index of se->nnE array from se->nnE*_index array based on neighbor type
		}

		neighbor_type = final_configuration[i];
		if (neighbor_type > 0) {
			++neighbor_type_cnt_final[neighbor_type - 1]; //number of type A/B/C neighboring atoms in final configuration
			//energy -= se->nnE[nnE_A_idxs[neighbor_type - 1]]*b_anisotropy_factor[i]; //A-A/B/C bond, told to leave in and comment out
		}
	}

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

	if ((neighbor_cnt_initial > 0) && (neighbor_cnt_final <= 1))
	{ // [ ]: this is preventing evaporation (when no final neighbors) 
		// BUG: this will mess up probabilities? by having fake transitions
		energy = 1000.;				// final configuration has no near neighbors,
									// so this effectively corresponds to an evaporation-like event.
									// Don't let it happen!
	}

	// ENHANCE: replace calculating the exp with memoizing up the value (uhash?) -> speedup?
	// overpotential is only for evaporation
	//*rate = 1e13*exp(-energy/(kBoltz*temperature)) //+ 1e-4*exp(-(energy-overpotential)/(kBoltz*temperature));
	*rate = 1e13*exp(-energy/(kBoltz*temperature));
	//printf("rate = %le\n", *rate);
	return 0;
}

/******************************************************************************/
/******************************************************************************/

int calculate_evaporation_rate(
	int initial_configuration[],			// initial configuration array of atom's nearest neighbors
	int atom_type,							// type of atom in consideration for transition
	double temperature,						// system temperature
	double overpotential,					// system overpotential
	double *rate,							// return value
	struct SimulationEnv* se
)
{
	int i;
	static double b_anisotropy_factor[12] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};			// optional anistropy factor
	double energy = 0.0;
	
	int neighbor_type; //type of atom for nearest neighbor
	int bond_idx, env_idx;
	if (se->solubility[atom_type])
	{
		//A can evaporate
		for (i=0;i<se->max_neighbors;++i)
		{ // calculate energy of initial state before evaporation
			neighbor_type = initial_configuration[i];
			bond_idx = get_bond_index(atom_type, neighbor_type, se);
			env_idx = get_env_index(1, bond_idx, se); // TODO: hard-coded 1st nn
			if (neighbor_type > 0)
				energy += se->nnEa[env_idx] * b_anisotropy_factor[i]; //bonding with A
		}
	}
	else
		// BUG: replace faking a really high energy with just not having the rate in the list
		energy = 1000.; //A cannot evaporate
	
	// dissolution/evaporation equation
	*rate = 1e4*exp(-(energy-overpotential)/(kBoltz*temperature));
		
	return 0;
}
