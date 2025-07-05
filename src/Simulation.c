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
	int transitioning_atom_idx, transition_jump_vector;
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
		//does this need to happen? - yes, overpotential/temperature changes
		
		// increment the elapsed time
		compute_transition_array(ss, se);
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
				k = ss->transition_probability.rate_arr_index[j];

				// pick the lucky atom
				// [ ]: third random number?
				// picks a type of transition and then which atom that has that transition will it act on?
				which_one = ss->rate_arr[k].transition_start_idx + (int)(drandj(&rand_seed) * (double)ss->rate_arr[k].transition_count);

				// which_one gives the location of the transition_arr, which gives
				// the info about the specific atom
				transitioning_atom_idx = ss->transition_arr[which_one]->atom_idx;
				transition_jump_vector = ss->transition_arr[which_one]->offset_idx;

				// if jump_vector == se->max_neighbors then the atom is going to evaporate
				moved_flag = true;

				adatom_before = 0;

				if (transition_jump_vector != se->max_neighbors)
				{
					//diffusion

					// coordinates atom is jumping to
					lastxt = ss->atom_arr[transitioning_atom_idx]->lattice[0] + jump_offset[transition_jump_vector].dx;
					lastyt = ss->atom_arr[transitioning_atom_idx]->lattice[1] + jump_offset[transition_jump_vector].dy;
					lastzt = ss->atom_arr[transitioning_atom_idx]->lattice[2] + jump_offset[transition_jump_vector].dz;

					adjust_pbc(&lastxt, &lastyt, &lastzt, se);

					atype = ss->atom_arr[transitioning_atom_idx]->type;

					// moves atom?
					remove_atom(transitioning_atom_idx, ss, se);
					natn = add_atom(lastxt, lastyt, lastzt, atype, NORMAL, ss, se);
				}
				else
				{
					// dissolution
					if (se->is_soluble[ss->atom_arr[transitioning_atom_idx]->type])
					{
						++ss->total_atoms_dissolved;
						remove_atom(transitioning_atom_idx, ss, se);	// evaporate the atom
						
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

		if (moved_flag == false) // ? only happens iff jump_vector == se->max_neighbors? (dissolution?)
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
void compute_transition_array(struct SimulationState* ss, struct SimulationEnv* se)
{
	
	int nonzero_rate_cnt = 0;
	double sum_of_rate_populations = 0.0;
	ss->frequency_sum = 0.0;
	
	int rate_const;
	Rate *r;
	for (int rate_idx=0; rate_idx < ss->rate_cnt; ++rate_idx)
	{
		r = &ss->rate_arr[rate_idx];
		if (r->transition_count != 0)
		{
			if (r->is_evaporation)
			{
				rate_const = calculate_evaporation_rate2(r->atom_env, ss->temperature, ss->overpotential, se);
			}
			else
			{
				rate_const = calculate_surf_diffusion_rate2(r->atom_env, ss->temperature, ss->overpotential, se);
			}
			r->k = rate_const;
			r->frequency = r->k*(double)r->transition_count;
			ss->frequency_sum += r->frequency;
			sum_of_rate_populations += r->transition_count;
			
			ss->transition_probability.rate_arr_index[nonzero_rate_cnt] = rate_idx;
			++nonzero_rate_cnt;
		}
	}

	// now compute bounds for jump probabilities
		
	int rate_idx;
	double current_probability = 0.0;

	for (int i=0; i<nonzero_rate_cnt; ++i)
	{
		ss->transition_probability.lbound[i] = current_probability;
		rate_idx = ss->transition_probability.rate_arr_index[i];
		current_probability += ss->rate_arr[rate_idx].frequency / ss->frequency_sum;
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
	
	unsigned char *atom_env = (unsigned char*) calloc(se->num_nn_types, sizeof(unsigned char));
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
		ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = -1;
   		next_x = ss->atom_arr[atom_idx]->lattice[0] + jump_offset[i].dx;
		next_y = ss->atom_arr[atom_idx]->lattice[1] + jump_offset[i].dy;
        next_z = ss->atom_arr[atom_idx]->lattice[2] + jump_offset[i].dz;

		adjust_pbc(&next_x, &next_y, &next_z, se);

        j = atom_at(next_x, next_y, next_z, ss->atom_arr, ss->zone_arr, se);

		// update atom_env
		if (j >= 0)
		{
			ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = j;
			
			int bond_idx = get_bond_index(ss->atom_arr[atom_idx]->type, ss->atom_arr[j]->type, se);
			int env_idx = get_env_index(1, bond_idx, se); // 1 for 1st nn
			
			atom_env[env_idx]++;
		}

		// // if no atom atom_idx that position but occ_neighbor array says there is, fix it
		// if ((j == -1) && (ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] >= 0))
		// 	ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = -1;
	}

	// cycle through the neighbor sites.  if there's an empty one, calculate the transition rate to it

	atom_rates_cnt = 0;

	// ENHANCE: redundant - already calculated initial config in j loop
	int intial_config_neighbor_cnt = get_initial_configuration(atom_idx, 0, se->max_neighbors, ss->atom_arr, start_config);		// k is number of near neighbors

	if (intial_config_neighbor_cnt == se->max_neighbors) // skip calculating a rate of a fully coordinated atom
		return atom_rates_cnt;

	for (int neighbor_idx=0; neighbor_idx < se->max_neighbors; ++neighbor_idx) // create transitions to each unoccupied neighbor
	{
		if (ss->atom_arr[atom_idx]->neighbor_atom_idxs[neighbor_idx] == -1) // if unoccupied, consider transition
		{ 
			// end_config is only used to identify the transitions that are functionally evaporations (no neighbors in end configuration)
			int final_config_neighbor_cnt = get_final_configuration(atom_idx, neighbor_idx, ss, se, end_config);

			if (final_config_neighbor_cnt < 1)
			{
				// transition in this direction leaves stray atom
				// these types are handled by evaporation
				continue;
			}
			
			// TODO: remove
			// calculate_surf_diffusion_rate(
			// 	start_config, // ENHANCE: make this look prettier
			// 	end_config,
			// 	ss->atom_arr[atom_idx]->type,
			// 	ss->temperature,
			// 	ss->overpotential,
			// 	&rate,
			// 	se
			// );
					
			++atom_rates_cnt;

			rate_idx = get_rate(atom_env, 0, ss, se);
			if (rate_idx == -1)
			{
				// if rate doesn't already exist, make new one
				rate_idx = create_new_rate(atom_env, 0, ss, se);
			}
			add_to_transition_list(rate_idx, atom_idx, neighbor_idx, ss, se);
		}
	}	

	// TODO: remove
	// calculate_evaporation_rate(	
	// 	start_config,
	// 	ss->atom_arr[atom_idx]->type,
	// 	ss->temperature,
	// 	ss->overpotential,
	// 	&rate,
	// 	se
	// );

	if (se->is_soluble[ss->atom_arr[atom_idx]->type])
	{
		// dissolution / evaporation transition
		rate_idx = get_rate(atom_env, 1, ss, se);
		if (rate_idx == -1)
		{
			// if rate doesn't already exist, make new one
			rate_idx = create_new_rate(atom_env, 1, ss, se);
		}
		// evaporation is considered to be last in jump_offset (not really in array but uses that index number)
		add_to_transition_list(rate_idx, atom_idx, se->max_neighbors, ss, se);
	}

	free(atom_env);
	atom_env = NULL;

	return atom_rates_cnt;						// gives number of current transitions for that atom
}

/******************************************************************************/
/******************************************************************************/

//is there a better way of doing this?
// checks if the rate constant is already in the rate list
int get_rate(unsigned char* atom_env, unsigned char is_evaporation, struct SimulationState* ss, struct SimulationEnv* se) 
{
	int env_cmp, evap_cmp;

	for (int i=0; i < ss->rate_cnt; ++i)
	{
		env_cmp = memcmp(ss->rate_arr[i].atom_env, atom_env, se->num_nn_types);
		env_cmp = env_cmp == 0;
		evap_cmp = ss->rate_arr[i].is_evaporation == is_evaporation;
		if (env_cmp && evap_cmp)
			return i;
	}

	return -1;
}

/******************************************************************************/
/******************************************************************************/
// create new Rate struct in rate array
// updates rate_array[rate_cnt], rate_cnt
int create_new_rate(unsigned char *atom_env, unsigned char is_evaporation, struct SimulationState* ss, struct SimulationEnv* se)
{
	// ss->rate_arr[ss->rate_cnt].k = rate;
	Rate *r = &ss->rate_arr[ss->rate_cnt];
	r->transition_start_idx = ss->transition_cnt;
	r->transition_count = 0;
	r->atom_env = (unsigned char *)malloc(se->num_nn_types * sizeof(unsigned char));
	memcpy(r->atom_env, atom_env, se->num_nn_types * sizeof(unsigned char));
	r->is_evaporation = is_evaporation;

	++ss->rate_cnt;

	return (ss->rate_cnt-1);
}

/******************************************************************************/
/******************************************************************************/

// add to rate_arr[rate_idx] the atom atom_idx going in direction offset_idx
// updates transition_arr, rate_arr[rate_idx].transition_count, atom_arr[atom_idx]->transition_indices[offset_idx]
void add_to_transition_list(int rate_idx, int atom_idx, int offset_idx, struct SimulationState* ss, struct SimulationEnv* se) // rate_arr index, atom_arr index, jump_offset index
{ // 0 2920 11
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
	{ // [ ]: what does this do? is this the same as in remove_transition? i=6
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
		ss->transition_arr[ss->transition_cnt] = NULL;
		// TODO: use pointers for rate array
		for (i = rate_idx + 1; i < ss->rate_cnt; ++i)
		{
			ss->rate_arr[i-1].transition_start_idx = ss->rate_arr[i].transition_start_idx;
			ss->rate_arr[i-1].transition_count = ss->rate_arr[i].transition_count;
			ss->rate_arr[i-1].k = ss->rate_arr[i].k;
			ss->rate_arr[i-1].frequency = ss->rate_arr[i].frequency;
			ss->rate_arr[i-1].atom_env = ss->rate_arr[i].atom_env;
			ss->rate_arr[i-1].is_evaporation = ss->rate_arr[i].is_evaporation;
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
	ss->transition_arr[ss->transition_cnt] = NULL;

	return;
}


/******************************************************************************/
/******************************************************************************/

void check_system(struct SimulationState* ss, struct SimulationEnv* se)
{
	int k, m, n, mm;
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

	for (int j = 0; j < ss->atom_cnt; ++j)
		for (int i=0;i<se->max_neighbors + 1;++i)				// extra 1 for evaporation
			{
				if (ss->atom_arr[j]->transition_indices[i] != -1)		// something can happen in the "i" direction
					take_off_transition_list(j, i, ss);
			}

	// for each atom, cycle through neighbor coordinates, and and reconcile occupancy

	for (int j = 0; j < ss->atom_cnt; ++j)
		for (int i=0; i < se->max_neighbors;++i)
		{
   			next_x = ss->atom_arr[j]->lattice[0] + jump_offset[i].dx;
			next_y = ss->atom_arr[j]->lattice[1] + jump_offset[i].dy;
	        next_z = ss->atom_arr[j]->lattice[2] + jump_offset[i].dz;

			adjust_pbc(&next_x, &next_y, &next_z, se);

			k = atom_at(next_x, next_y, next_z, ss->atom_arr, ss->zone_arr, se);

			if (k >= 0)
			{
				// an atom has been found atom_idx this neighbor site.

				ss->atom_arr[j]->neighbor_atom_idxs[i] = k;
				ss->atom_arr[k]->neighbor_atom_idxs[opposite_offset[i]] = j;
			}
		}

	// now we'll reconcile buried atoms.
	// does this process need to happen if nothing is buried???
	do
	{
		errors = 0;

		for (int j = 0; j < ss->atom_cnt; ++j)
		for (int i=0;i<se->max_neighbors;++i)
		{
			if (ss->atom_arr[j]->neighbor_atom_idxs[i] == -2)
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

						n = ss->atom_arr[m]->neighbor_atom_idxs[opposite_offset[k]];
						if (n != -2)
						{
							if (n == -1) 
								ss->atom_arr[m]->neighbor_atom_idxs[opposite_offset[k]] = -2;
							if (n == -3)
								ss->atom_arr[j]->neighbor_atom_idxs[i] = -3;		// random trumps buried

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

		for (int j = 0; j < ss->atom_cnt; ++j)
			for (int i = 0; i < se->max_neighbors; ++i)
			{
				if (ss->atom_arr[j]->neighbor_atom_idxs[i] == -3)
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

							n = ss->atom_arr[m]->neighbor_atom_idxs[opposite_offset[k]];
							if (n != -3)
							{
								ss->atom_arr[m]->neighbor_atom_idxs[opposite_offset[k]] = -3;
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
			if (atom[j]->neighbor_atom_idxs[i] >= 0)
			{	
				if (atom[j]->type == atom[atom[j]->neighbor_atom_idxs[i]]->type)
				++k;
			}

			if (atom[j]->neighbor_atom_idxs[i] == -2)
				++k;
		}

		// SPECIAL SC routine included here otherwise for second nearest neighbors?

		if (k == se->max_neighbors)
		{
			// bury atom j

			for (i=0;i<se->max_neighbors;++i)
			{
				m = atom[j]->neighbor_atom_idxs[i];
				if (m >= 0)
					atom[m]->neighbor_atom_idxs[opposite_offset[i]] = -2;
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
			atom[atom_cnt-1] = NULL;
			--atom_cnt;
			--j;
		}
	}*/

	// now we can recalculate diffusion rates

	for (int j = 0; j < ss->atom_cnt; ++j)
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
		neighbor_type = initial_configuration[i];

		bond_idx = get_bond_index(atom_type, neighbor_type, se);
		env_idx = get_env_index(1, bond_idx, se); // TODO: hard-coded 1st nn

		if (neighbor_type >= 0) {
			++neighbor_type_cnt_initial[neighbor_type]; //number of type A/B/C neighboring atoms in init configuration
			energy += se->nn_energy[env_idx] * b_anisotropy_factor[i]; //A-A/B/C bond - grab correct index of se->nnE array from se->nnE*_index array based on neighbor type
		}

		neighbor_type = final_configuration[i];
		if (neighbor_type >= 0) {
			++neighbor_type_cnt_final[neighbor_type]; //number of type A/B/C neighboring atoms in final configuration
			//energy -= se->nnE[nnE_A_idxs[neighbor_type]]*b_anisotropy_factor[i]; //A-A/B/C bond, told to leave in and comment out
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

// calculates surface diffusion rate using atom environment
double calculate_surf_diffusion_rate2(
	unsigned char* atom_env,
	double temperature,						// system temperature
	double overpotential,					// system overpotential
	struct SimulationEnv* se
)
{ // ENHANCE: pass the number of nearest neighbors? (best practice)
	int i;
	double energy = 0.0;
	
	// why is this static, not const? it should also be a simulation input
	// "optional" anistropy factor, indices match jump_offset
	// static double b_anisotropy_factor[12] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; 

	int neighbor_cnt_initial = 0;
	for (int i = 0; i < se->num_nn_types; i++)
	{
		energy += se->nn_energy[i] * atom_env[i];
		// TODO: anisotropy factor would require storing directions in Rate (probably as sum, to be multiplied with energy)
		neighbor_cnt_initial++;
	}

	// these override the previous energy sum
	if (neighbor_cnt_initial == 0)
	{
		// no neighbors - this condition corresponds to a diffuser walking through a lattice (a lattice gas)
		energy = -1.0;
	}

	// ENHANCE: replace calculating the exp with memoizing up the value (uhash?) -> speedup?
	// overpotential is only for evaporation
	//*rate = 1e13*exp(-energy/(kBoltz*temperature)) //+ 1e-4*exp(-(energy-overpotential)/(kBoltz*temperature));
	return 1e13*exp(-energy/(kBoltz*temperature));
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
	if (se->is_soluble[atom_type])
	{
		//A can evaporate
		for (i=0;i<se->max_neighbors;++i)
		{ // calculate energy of initial state before evaporation
			neighbor_type = initial_configuration[i];
			bond_idx = get_bond_index(atom_type, neighbor_type, se);
			env_idx = get_env_index(1, bond_idx, se); // TODO: hard-coded 1st nn
			if (neighbor_type >= 0)
				energy += se->nn_energy[env_idx] * b_anisotropy_factor[i]; //bonding with A
		}
	}
	else
		// BUG: replace faking a really high energy with just not having the rate in the list
		energy = 1000.; //A cannot evaporate
	
	// dissolution/evaporation equation
	*rate = 1e4*exp(-(energy-overpotential)/(kBoltz*temperature));
		
	return 0;
}

// calculates evaporation rate using atom environment
double calculate_evaporation_rate2(
	unsigned char *atom_env,
	double temperature,						// system temperature
	double overpotential,					// system overpotential
	struct SimulationEnv* se
)
{
	double energy = 0.0;
	
	static double b_anisotropy_factor[12] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};			// optional anistropy factor

	for (int i = 0; i < se->num_nn_types; i++)
	{
		energy += se->nn_energy[i] * atom_env[i];
		// TODO: anisotropy factor would require storing directions in Rate (probably as sum, to be multiplied with energy)
	}
	
	// dissolution/evaporation equation
	return 1e4*exp(-(energy-overpotential)/(kBoltz*temperature));
		
	return 0;
}
