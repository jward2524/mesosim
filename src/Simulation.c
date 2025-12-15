#include "Simulation.h"
#include "Random.h"
#include "Simulation_Aux.h"
#include "FileIO.h"
#include "Atoms.h"
#include "Mesosim.h"
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

int rate_skip; // transition list binary search index
int adatom_before; // XXX: never used?

int lastxt, lastyt, lastzt; // containers for coordinates of a next step
bool checkpoint_reached = false;

const double FABS_TOL = 1e-6;

// ENHANCE: pass struct with all simulation parameters as argument
unsigned long perform_simulation(struct SimulationState* ss, struct SimulationEnv* se, struct LoggingState* ls) //potentially FILE* as arguments
{
	//printf("Simulation starting\n");
	// long int i;
	int j, k;

	double cur_stime, prev_stime; // prev and current times, for overpotential moving

	double transition_type_probability;

	//double vap;

	int selected_transition_idx;
	int transitioning_atom_idx, transition_jump_vector;
	int transitioned_atom_idx; //changed from nan bc keyword
	int atype;

	int moved_flag = true;

	// ls->framenum = 0; 

	prev_stime = ss->elapsed_stime;
	
	// ss->total_atoms_dissolved = 0;

	// // iteration count
	// ss->iter = 0;
	
	bool simulation_end = false;
	
	//writes data time intervals and run time in original code

	//simulation_is_going = true;

	//print system and zone sizes to the file?

	/*if (num_sims > 0) //this probably will not happen
	{
		initialize_initial_structure(simulation_type);
		elapsed_stime = 0.0;
	}*/ 

	// initialize simulation kinetics, and draw a picture

	//printf("transition time!\n");
	for (int i=0; i < ss->atom_cnt; ++i)
	{
		// resets all kinetic paramters
		refresh_transitions(i, ss, se);
	}
	compute_transition_array(ss, se);
	
	// TODO: move this out from here, to only where output is necessary
	organize(ss->atom_arr, ss->atom_cnt, se->primitive_basis); //replacement for copy_xyz_to_coord but might not be necessary

	char suffix[BUFFER_SIZE];
	if (ls->analysis_type == REGULAR_TIME_INTERVALS)
	{
		snprintf(suffix, BUFFER_SIZE, "t0");
	}
	else if (ls->analysis_type == LN_TIME_INTERVALS)
	{
		snprintf(suffix, BUFFER_SIZE, "t0");
	}
	else if (ls->analysis_type == ITERATION_INTERVALS)
	{
		snprintf(suffix, BUFFER_SIZE, "i0");
	}

	// initial state
	output_log_file(ls->sim_log_file, ls->framenum, ss->iter, ss->elapsed_stime, ss->temperature, ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
	write_xyz_file(ls->position_log_prefix, ls->framenum, suffix, ss, se);
	ls->framenum++;

	// start with everything current to the current state
	// choose a transition
	// perform the transition
	// update time
	// update rates
	// next iteration

	int old_x, old_y, old_z;
	// int neighbor_x, neighbor_y, neighbor_z, neighbor_idx;

	printf("Setup complete, iteration start");
	while (!simulation_end)
	{
		ss->iter++;

		if (ss->simulation_should_kill_itself) // abort simulation (only happens if atoms overlap)
		{
			output_log_file(ls->sim_log_file, ls->framenum, ss->iter, ss->elapsed_stime, ss->temperature, ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
			write_xyz_file(ls->position_log_prefix, ls->framenum, suffix, ss, se);

			ss->simulation_should_kill_itself = false;

			// TODO: remove from here
			// organize(ss->atom_arr, ss->atom_cnt, se->primitive_basis); // XXX: remove from here
			printf("killed somehow\n");

			return 1; //return 1 b/c error
		}

		// pick the type of transition to occur
		transition_type_probability = drandj(&rand_seed);
		rate_skip = ss->rate_cnt / 2;
		j = rate_skip;

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
				selected_transition_idx = ss->rate_arr[k].transition_start_idx + (int)(drandj(&rand_seed) * (double)ss->rate_arr[k].transition_count);

				// selected_transition_idx gives the location of the transition_arr, which gives
				// the info about the specific atom
				transitioning_atom_idx = ss->transition_arr[selected_transition_idx]->atom_idx;
				transition_jump_vector = ss->transition_arr[selected_transition_idx]->offset_idx;

				// if jump_vector == se->num_transition_vectors then the atom is going to evaporate
				moved_flag = true;

				adatom_before = 0;

				old_x = ss->atom_arr[transitioning_atom_idx]->lattice[0];
				old_y = ss->atom_arr[transitioning_atom_idx]->lattice[1];
				old_z = ss->atom_arr[transitioning_atom_idx]->lattice[2];

				// perform transition
				if (transition_jump_vector != se->num_transition_vectors)
				{
					//diffusion

					// coordinates atom is jumping to
					lastxt = ss->atom_arr[transitioning_atom_idx]->lattice[0] + se->transition_vectors[transition_jump_vector].dx;
					lastyt = ss->atom_arr[transitioning_atom_idx]->lattice[1] + se->transition_vectors[transition_jump_vector].dy;
					lastzt = ss->atom_arr[transitioning_atom_idx]->lattice[2] + se->transition_vectors[transition_jump_vector].dz;

					adjust_pbc(&lastxt, &lastyt, &lastzt, se);

					atype = ss->atom_arr[transitioning_atom_idx]->type;

					// moves atom?
					remove_atom(transitioning_atom_idx, ss, se);
					transitioned_atom_idx = add_atom(lastxt, lastyt, lastzt, atype, NORMAL, ss, se);
				}
				else
				{
					// dissolution

					// if not soluble, then there was an issue somewhere
					if (!(se->is_soluble[ss->atom_arr[transitioning_atom_idx]->type]))
					{
						fprintf(stderr, "Attempting to dissolve an insoluble atom - Terminatin\n");
       					clean_and_exit(errno);
					}
					++ss->total_atoms_dissolved;
					remove_atom(transitioning_atom_idx, ss, se);	// evaporate the atom
					transitioned_atom_idx = -1;
				}
				// end of performing transition
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
			else
			{
				fprintf(stderr, "Transition decision failed\n");
				clean_and_exit(1);
			}
		}

		if (moved_flag == false) // ? only happens iff jump_vector == se->num_transition_vectors? (dissolution?)
		{
			fprintf(stderr, "for some reason I didn't transition\n");
			clean_and_exit(1);
		}
		
		// increment the elapsed time (how much passed before last transition)
		// time passed, then transition occurred
		// whether this is the 1st, 2nd, or 3rd random number in the loop affects the time (diff random numbers)
		ss->elapsed_stime -= log(drandj(&rand_seed)) / ss->frequency_sum;
		
		if (se->overpotential_ramp_rate != 0.0)
		{
			cur_stime = ss->elapsed_stime;
			ss->overpotential += (cur_stime-prev_stime) * se->overpotential_ramp_rate;
			prev_stime = ss->elapsed_stime;
		}
		
		// update rates
		update_outdated_transitions(old_x, old_y, old_z, transitioned_atom_idx, ss, se);
		
		// ENHANCE - not all are necessary if overpotential/temperature don't change?
		compute_transition_array(ss, se);
		// end of updating rates
	
		// after iteration, log if necessary
		if (ss->iter % 100 == 0)
			printf("iteration %ld, time %le\n", ss->iter, ss->elapsed_stime);

		// TODO: implement the checkpoint lists
		if ((ls->analysis_type == REGULAR_TIME_INTERVALS) || (ls->analysis_type == LN_TIME_INTERVALS))\
		{
			checkpoint_reached = (ss->elapsed_stime >= ls->next_log_checkpoint);
			snprintf(suffix, BUFFER_SIZE, "0");
		}
		else if (ls->analysis_type == ITERATION_INTERVALS)
		{
			checkpoint_reached = fabs(ls->next_log_checkpoint - (double) ss->iter) < FABS_TOL;
			if (!checkpoint_reached && (ss->iter > ls->next_log_checkpoint))
			{
				fprintf(stderr, "Iterations (%lu) exceeded log checkpoint (%lf) without noticing", ss->iter, ls->next_log_checkpoint);
				clean_and_exit(1);
			}
		}

		if (checkpoint_reached)
		{
			// organize(ss->atom_arr, ss->atom_cnt, se->primitive_basis); //replaced but do i really need it
			
			// update suffix and increment next checkpoint 
			if (ls->analysis_type == REGULAR_TIME_INTERVALS)
			{
				snprintf(suffix, BUFFER_SIZE, "t%.4lf", ls->next_log_checkpoint);
				
				// bring next_log_checkpoint up to and one step beyond elapsed_stime
				while (ls->next_log_checkpoint <= ss->elapsed_stime)
				{
					ls->next_log_checkpoint += ls->log_interval;
				}
			}
			else if (ls->analysis_type == LN_TIME_INTERVALS)
			{
				snprintf(suffix, BUFFER_SIZE, "t%.4lf", ls->next_log_checkpoint);
				while (ls->next_log_checkpoint <= ss->elapsed_stime)
				{
					ls->next_log_checkpoint *= ls->log_interval;
				}
			}
			else if (ls->analysis_type == ITERATION_INTERVALS)
			{
				snprintf(suffix, BUFFER_SIZE, "i%lu", (unsigned long) ls->next_log_checkpoint);
				ls->next_log_checkpoint += ls->log_interval;
			}

			// record the elapsed time in a file here
			printf("writing file %d: elapsed_stime = %le\n", ls->framenum, ss->elapsed_stime);
			output_log_file(ls->sim_log_file, ls->framenum, ss->iter, ss->elapsed_stime, ss->temperature, ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
			write_xyz_file(ls->position_log_prefix, ls->framenum, suffix, ss, se);

			++ls->framenum;
		}

		// check if simulation is over
		if (ss->sim_end_type == SIM_END_BY_STIME){
			simulation_end = (ss->elapsed_stime >= ss->run_stime);
		}
		else if (ss->sim_end_type == SIM_END_BY_ITERATIONS) {
			simulation_end = (ss->iter >= ss->final_iteration);
		}
	}

	//write elapsed_stime to mark finish
	output_log_file(ls->sim_log_file, ls->framenum, ss->iter, ss->elapsed_stime, ss->temperature, ss->overpotential, ss->atom_cnt, ss->total_internal_energy);
	
	snprintf(suffix, BUFFER_SIZE, "%s_final", suffix);
	write_xyz_file(ls->position_log_prefix, ls->framenum, suffix, ss, se);
	
	if ((ss->final_iteration > 0) && (ss->iter >= ss->final_iteration))
	{
		fprintf(ls->sim_log_file, "Reached final iteration and terminated\n");
	}
	if ((ss->run_stime > 0) && (ss->elapsed_stime >= ss->run_stime))
	{
		fprintf(ls->sim_log_file, "Reached end of simulation time and terminated\n");
	}

	printf("Finished simulation\n");
		
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
	
	double rate_const;
	Rate *r;
	// ENHANCE: parallelize
	for (int rate_idx=0; rate_idx < ss->rate_cnt; ++rate_idx)
	{
		r = &ss->rate_arr[rate_idx];
		if (r->transition_count != 0)
		{
			if (r->is_evaporation)
			{
				rate_const = calculate_evaporation_rate(r->atom_env, ss->temperature, ss->overpotential, se);
			}
			else
			{
				rate_const = calculate_surf_diffusion_rate(r->atom_env, r->final_config_neighbor_cnt, ss->temperature, se);
			}
			if (isnan(rate_const))
			{
				fprintf(stderr, "Rate constant is nan");
				clean_and_exit(1);
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

void update_outdated_transitions(int old_x, int old_y, int old_z, int transitioned_atom_idx, struct SimulationState* ss, struct SimulationEnv* se)
{
	// transitioned_atom_index - new atom index (-1 if evaporated)
	int neighbor_x, neighbor_y, neighbor_z, neighbor_idx;
	
	// refresh neighbors of previous site
	for (int i = 0; i < se->num_energy_contributors; i++)
	{
		neighbor_x = old_x + se->transition_vectors[i].dx;
		neighbor_y = old_y + se->transition_vectors[i].dy;
		neighbor_z = old_z + se->transition_vectors[i].dz;

		adjust_pbc(&neighbor_x, &neighbor_y, &neighbor_z, se);

		neighbor_idx = atom_at(neighbor_x, neighbor_y, neighbor_z, ss->atom_arr, ss->zone_arr, se);

		if (neighbor_idx >= 0)
			refresh_transitions(neighbor_idx, ss, se);
	}

	// if transition was not evaporation
	if (transitioned_atom_idx >= 0)
	{
		// refresh moved atom - also gets refreshed as neighbor of previous site?
		refresh_transitions(transitioned_atom_idx, ss, se);

		// refresh neighbors of moved-to site
		for (int i = 0; i < se->num_energy_contributors; i++)
		{
			neighbor_idx = ss->atom_arr[transitioned_atom_idx]->neighbor_atom_idxs[i];

			if (neighbor_idx >= 0)
				refresh_transitions(neighbor_idx, ss, se);
		}
	}
}

/******************************************************************************/
/******************************************************************************/
// updates [atom_arr[atom_idx], neighbor_atom_idxs, transition_arr[i], rate_arr[i].transition_start_idx], initializes rate_arr
int refresh_transitions(int atom_idx, struct SimulationState* ss, struct SimulationEnv* se) // atom_idx = index on atom list
{
	int neighbor_idx;
	int rate_idx; // position of rate rate in rate list rate_arr[]
	int next_x, next_y, next_z;

	int atom_rates_cnt;	// this is returned as the number of transitions (se->jump_offsets) this atom can undergo, excluding evaporation
	int start_config[MAXIMUM_NUMBER_OF_NEIGHBORS]; // -1 if empty, type if filled
	int end_config[MAXIMUM_NUMBER_OF_NEIGHBORS];
	
	unsigned char *atom_env = (unsigned char*) calloc(se->num_nn_types, sizeof(unsigned char));
	// unsigned char env_hash[se->num_nn_levels * se->num_bond_types];

	// first, remove all mention of this atom from transition list
	// extra 1 for evaporation
	for (int i=0; i < se->num_transition_vectors + se->dissolution; ++i) 
	{
		// transition can happen in the "i" direction
		if (ss->atom_arr[atom_idx]->transition_indices[i] != -1) 
			take_off_transition_list(atom_idx, i, ss);
	}

	// update neighbors
	// cycle through neighbor coordinates, and check if there is an atom there
	for (int i=0; i < se->num_transition_vectors; ++i)
	{
		ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = -1;

   		next_x = ss->atom_arr[atom_idx]->lattice[0] + se->transition_vectors[i].dx;
		next_y = ss->atom_arr[atom_idx]->lattice[1] + se->transition_vectors[i].dy;
        next_z = ss->atom_arr[atom_idx]->lattice[2] + se->transition_vectors[i].dz;

		adjust_pbc(&next_x, &next_y, &next_z, se);

        neighbor_idx = atom_at(next_x, next_y, next_z, ss->atom_arr, ss->zone_arr, se);

		if (neighbor_idx >= 0)
		{
			ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = neighbor_idx;
		}

		// // if no atom atom_idx that position but occ_neighbor array says there is, fix it
		// if ((j == -1) && (ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] >= 0))
		// 	ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] = -1;
	}

	// update atom_env and energy
	ss->total_internal_energy -= ss->atom_arr[atom_idx]->energy;
	ss->atom_arr[atom_idx]->energy = 0;

	for (int i=0; i < se->num_energy_contributors; ++i)
	{
   		next_x = ss->atom_arr[atom_idx]->lattice[0] + se->transition_vectors[i].dx;
		next_y = ss->atom_arr[atom_idx]->lattice[1] + se->transition_vectors[i].dy;
        next_z = ss->atom_arr[atom_idx]->lattice[2] + se->transition_vectors[i].dz;

		adjust_pbc(&next_x, &next_y, &next_z, se);

        neighbor_idx = atom_at(next_x, next_y, next_z, ss->atom_arr, ss->zone_arr, se);

		// update atom_env
		if (neighbor_idx >= 0)
		{
			int nn_level;
			int diff = i;
			for (int j=0; j<se->num_nn_levels; j++)
			{
				diff = diff - se->atoms_per_nn_level[j];
				if (diff < 0)
				{
					nn_level = j;
					break;
				}
			}
			int bond_idx = get_bond_index(ss->atom_arr[atom_idx]->type, ss->atom_arr[neighbor_idx]->type, se);
			int env_idx = get_env_index(nn_level, bond_idx, se);
			
			atom_env[env_idx]++;
			ss->atom_arr[atom_idx]->energy += se->nn_energy[env_idx];
		}
	}

	ss->atom_arr[atom_idx]->energy /= 2;
	ss->total_internal_energy += ss->atom_arr[atom_idx]->energy;

	// cycle through the neighbor sites.  if there's an empty one, calculate the transition rate to it
	atom_rates_cnt = 0;

	// ENHANCE: redundant - already calculated initial config in j loop
	int intial_config_neighbor_cnt = get_initial_configuration(atom_idx, se->num_transition_vectors, ss->atom_arr, start_config);		// k is number of near neighbors

	// skip calculating a rate of a fully coordinated atom
	if (intial_config_neighbor_cnt == se->num_transition_vectors) 
		return atom_rates_cnt;

	// create transitions to each unoccupied neighbor
	int final_config_neighbor_cnt;
	for (int i=0; i < se->num_transition_vectors; ++i) 
	{
		// if unoccupied, consider transition
		if (ss->atom_arr[atom_idx]->neighbor_atom_idxs[i] == -1) 
		{ 
			// end_config is only used to identify the transitions that are functionally evaporations (no neighbors in end configuration)
			final_config_neighbor_cnt = get_final_configuration(atom_idx, i, ss, se, end_config);

			if (final_config_neighbor_cnt <= 1)
			{
				// transition in this direction leaves stray atom
				// or leads to pair climbs
				// these types are handled by evaporation
				continue;
			}

			++atom_rates_cnt;

			rate_idx = get_rate(atom_env, final_config_neighbor_cnt, 0, ss, se);
			if (rate_idx == -1)
			{
				// if rate doesn't already exist, make new one
				rate_idx = create_new_rate(atom_env, final_config_neighbor_cnt, 0, ss, se);
			}
			add_to_transition_list(rate_idx, atom_idx, i, ss, se);
		}
	}	

	// if atom type is soluble, add dissolution transition
	if (se->is_soluble[ss->atom_arr[atom_idx]->type])
	{
		// dissolution / evaporation transition
		final_config_neighbor_cnt = -1;
		rate_idx = get_rate(atom_env, final_config_neighbor_cnt, 1, ss, se);
		if (rate_idx == -1)
		{
			// if rate doesn't already exist, make new one
			rate_idx = create_new_rate(atom_env, final_config_neighbor_cnt, 1, ss, se);
		}
		// evaporation is considered to be last in se->transition_vectors (not really in array but uses that index number)
		add_to_transition_list(rate_idx, atom_idx, se->num_transition_vectors, ss, se);
	}

	free(atom_env);
	atom_env = NULL;

	return atom_rates_cnt;						// gives number of current transitions for that atom
}

/******************************************************************************/
/******************************************************************************/

//is there a better way of doing this?
// checks if the rate constant is already in the rate list
int get_rate(unsigned char* atom_env, int final_config_neighbor_cnt, unsigned char is_evaporation, struct SimulationState* ss, struct SimulationEnv* se) 
{
	int env_cmp, evap_cmp, finalcnt_cmp;

	for (int i=0; i < ss->rate_cnt; ++i)
	{
		env_cmp = memcmp(ss->rate_arr[i].atom_env, atom_env, se->num_nn_types);
		env_cmp = env_cmp == 0;
		evap_cmp = ss->rate_arr[i].is_evaporation == is_evaporation;
		finalcnt_cmp = ss->rate_arr[i].final_config_neighbor_cnt == final_config_neighbor_cnt;
		if (env_cmp && evap_cmp && finalcnt_cmp)
			return i;
	}

	return -1;
}

/******************************************************************************/
/******************************************************************************/
// create new Rate struct in rate array
// updates rate_array[rate_cnt], rate_cnt
int create_new_rate(unsigned char *atom_env, int final_config_neighbor_cnt, unsigned char is_evaporation, struct SimulationState* ss, struct SimulationEnv* se)
{
	// ss->rate_arr[ss->rate_cnt].k = rate;
	Rate *r = &ss->rate_arr[ss->rate_cnt];
	r->transition_start_idx = ss->transition_cnt;
	r->transition_count = 0;
	r->atom_env = (unsigned char *)malloc(se->num_nn_types * sizeof(unsigned char));
	memcpy(r->atom_env, atom_env, se->num_nn_types * sizeof(unsigned char));
	r->is_evaporation = is_evaporation;

	++ss->rate_cnt;
	if (ss->rate_cnt > se->max_rates)
	{
		fprintf(stderr, "Number of rates (%lld) is exceeding set maximum (%lld)", ss->rate_cnt, se->max_rates);
        clean_and_exit(errno);
	}

	return (ss->rate_cnt-1);
}

/******************************************************************************/
/******************************************************************************/

// add to rate_arr[rate_idx] the atom atom_idx going in direction offset_idx
// updates transition_arr, rate_arr[rate_idx].transition_count, atom_arr[atom_idx]->transition_indices[offset_idx]
void add_to_transition_list(int rate_idx, int atom_idx, int offset_idx, struct SimulationState* ss, struct SimulationEnv* se) // rate_arr index, atom_arr index, se->transition_vectors index
{ // 0 2920 11
	int i; // loop variable
	int n;
	int initial_transition_index, final_transition_index; // initial and final transition_arr index

	// make room for the new arrival
	ss->transition_arr[ss->transition_cnt] = (Transition *)malloc(sizeof(Transition));	// adds entry to the end of the list
	if (ss->transition_arr[ss->transition_cnt] == NULL)
	{
		// TODO: free mallocs before exiting
		fprintf(stderr, "Couldn't allocate memory for atom %lld: %s\n", ss->transition_cnt, strerror(errno));
		clean_and_exit(errno);
	}
	if ((unsigned int) ss->transition_cnt > se->max_transitions)
	{
		fprintf(stderr, "More transitions (%lld) than allocated in transition array (%llu)\n", ss->transition_cnt, se->max_transitions);
		clean_and_exit(errno);
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
	if (ss->transition_cnt > se->max_transitions)
	{
		fprintf(stderr, "Number of transitions (%lld) is exceeding set maximum (%lld)", ss->transition_cnt, se->max_transitions);
        clean_and_exit(errno);
	}

	return;
}

/******************************************************************************/
/******************************************************************************/
// updates atom_arr[atom_idx], transition_arr[i], rate_arr[i].transition_start_idx
void take_off_transition_list(int atom_idx, int offset_idx, struct SimulationState* ss)	// removes atom jumping in the se->transition_vectors[offset_idx] direction
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
	if (ss->atom_cnt < 0)
	{
		fprintf(stderr, "Number of transitions (%lld) has dropped below zero", ss->transition_cnt);
        clean_and_exit(errno);
	}

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
		if (ss->rate_cnt < 0)
		{
			fprintf(stderr, "Number of rates (%lld) has dropped below zero", ss->rate_cnt);
			clean_and_exit(errno);
		}

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
	int k, m, n;
	int errors;

	int next_x, next_y, next_z;

	int nnx, nny, nnz;

	// does a careful check to make sure that a system is ready to be simulated.
	// assumptions:  (1) all real atoms are in the places they think they are
	// (2) all buried atoms are actually buried.

	// first, remove all atoms from the transition list.  We'll add them after we check neighbors

	for (int j = 0; j < ss->atom_cnt; ++j)
		for (int i=0;i<se->num_transition_vectors + 1;++i)				// extra 1 for evaporation
			{
				if (ss->atom_arr[j]->transition_indices[i] != -1)		// something can happen in the "i" direction
					take_off_transition_list(j, i, ss);
			}

	// for each atom, cycle through neighbor coordinates, and and reconcile occupancy

	for (int j = 0; j < ss->atom_cnt; ++j)
		for (int i=0; i < se->num_transition_vectors;++i)
		{
   			next_x = ss->atom_arr[j]->lattice[0] + se->transition_vectors[i].dx;
			next_y = ss->atom_arr[j]->lattice[1] + se->transition_vectors[i].dy;
	        next_z = ss->atom_arr[j]->lattice[2] + se->transition_vectors[i].dz;

			adjust_pbc(&next_x, &next_y, &next_z, se);

			k = atom_at(next_x, next_y, next_z, ss->atom_arr, ss->zone_arr, se);

			if (k >= 0)
			{
				// an atom has been found atom_idx this neighbor site.

				ss->atom_arr[j]->neighbor_atom_idxs[i] = k;
				ss->atom_arr[k]->neighbor_atom_idxs[se->opposite_tvectors[i]] = j;
			}
		}

	// now we'll reconcile buried atoms.
	// does this process need to happen if nothing is buried???
	do
	{
		errors = 0;

		for (int j = 0; j < ss->atom_cnt; ++j)
		for (int i=0;i<se->num_transition_vectors;++i)
		{
			if (ss->atom_arr[j]->neighbor_atom_idxs[i] == -2)
			{
				// find coordinate of buried atom

				next_x = ss->atom_arr[j]->lattice[0] + se->transition_vectors[i].dx;
				next_y = ss->atom_arr[j]->lattice[1] + se->transition_vectors[i].dy;
		        next_z = ss->atom_arr[j]->lattice[2] + se->transition_vectors[i].dz;

				adjust_pbc(&next_x, &next_y, &next_z, se);

				for (k=0;k<se->num_transition_vectors;++k)
				{
					nnx = next_x + se->transition_vectors[k].dx;
					nny = next_y + se->transition_vectors[k].dy;
					nnz = next_z + se->transition_vectors[k].dz;

					adjust_pbc(&nnx, &nny, &nnz, se);

					m = atom_at(nnx, nny, nnz, ss->atom_arr, ss->zone_arr, se);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = ss->atom_arr[m]->neighbor_atom_idxs[se->opposite_tvectors[k]];
						if (n != -2)
						{
							if (n == -1) 
								ss->atom_arr[m]->neighbor_atom_idxs[se->opposite_tvectors[k]] = -2;
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
			for (int i = 0; i < se->num_transition_vectors; ++i)
			{
				if (ss->atom_arr[j]->neighbor_atom_idxs[i] == -3)
				{
					// find coordinate of buried atom

					next_x = ss->atom_arr[j]->lattice[0] + se->transition_vectors[i].dx;
					next_y = ss->atom_arr[j]->lattice[1] + se->transition_vectors[i].dy;
					next_z = ss->atom_arr[j]->lattice[2] + se->transition_vectors[i].dz;

					adjust_pbc(&next_x, &next_y, &next_z, se);

					for (k=0;k<se->num_transition_vectors;++k)
					{
						nnx = next_x + se->transition_vectors[k].dx;
						nny = next_y + se->transition_vectors[k].dy;
						nnz = next_z + se->transition_vectors[k].dz;

						adjust_pbc(&nnx, &nny, &nnz, se);

						m = atom_at(nnx, nny, nnz, ss->atom_arr, ss->zone_arr, se);

						if ((m >= 0)&&(m!= j))
						{
							// another atom (m) is connected to this atom.  If it sees this position as 
							// a buried atom, great.  Otherwise, reconcile

							n = ss->atom_arr[m]->neighbor_atom_idxs[se->opposite_tvectors[k]];
							if (n != -3)
							{
								ss->atom_arr[m]->neighbor_atom_idxs[se->opposite_tvectors[k]] = -3;
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

		for (i=0;i<se->num_transition_vectors;++i)
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

		if (k == se->num_transition_vectors)
		{
			// bury atom j

			for (i=0;i<se->num_transition_vectors;++i)
			{
				m = atom[j]->neighbor_atom_idxs[i];
				if (m >= 0)
					atom[m]->neighbor_atom_idxs[se->opposite_tvectors[i]] = -2;
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
// calculates surface diffusion rate using atom environment
double calculate_surf_diffusion_rate(
	unsigned char* atom_env,
	int final_config_neighbor_cnt,
	double temperature,						// system temperature
	// double overpotential,					// system overpotential
	struct SimulationEnv* se
)
{ 
	// ENHANCE: pass the number of nearest neighbors? (best practice)
	double energy = 0.0;
	
	// why is this static, not const? it should also be a simulation input
	// "optional" anistropy factor, indices match se->transition_vectors
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
	
	// to simulate Ehrlich-Schwoebel barrier
	if (final_config_neighbor_cnt <= 2)
	{
		energy = 100;
	}

	// ENHANCE: replace calculating the exp with memoizing up the value (uhash?) -> speedup?
	// overpotential is only for evaporation
	//*rate = 1e13*exp(-energy/(kBoltz*temperature)) //+ 1e-4*exp(-(energy-overpotential)/(kBoltz*temperature));
	return 1e13*exp(-energy/(kBoltz*temperature));
}

/******************************************************************************/
/******************************************************************************/
// calculates evaporation rate using atom environment
double calculate_evaporation_rate(
	unsigned char *atom_env,
	double temperature,						// system temperature
	double overpotential,					// system overpotential
	struct SimulationEnv* se
)
{
	double energy = 0.0;
	
	// static double b_anisotropy_factor[12] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};			// optional anistropy factor

	for (int i = 0; i < se->num_nn_types; i++)
	{
		energy += se->nn_energy[i] * atom_env[i];
		// TODO: anisotropy factor would require storing directions in Rate (probably as sum, to be multiplied with energy)
	}
	
	// dissolution/evaporation equation
	return 1e4*exp(-(energy-overpotential)/(kBoltz*temperature));
		
	return 0;
}
