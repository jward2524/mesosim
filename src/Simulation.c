#include "stdafx.h"
#include "Defs.h"
#include "Geometry.h"
#include "Vector.h"
#include "Random.h"
#include "Simulation.h"
#include "Simulation_Aux.h"
#include "Atoms.h"
#include "FileIO.h"

int rate_skip;

double xr, yr, zr;
int rwx, rwy, rwz;
double rw[3], rrp[3];
double ta1, ta2;

int adatom_before;

char atom_names[3][3]={"1", "2", "3"};
double default_color[3] = {0., 0., 0.};

int simulation_type = SIMULATION_TYPE_UNDEFINED;
int atom_cnt = 0;
Atom temp_atom;
bool simulation_should_kill_itself;
Atom* atom_arr[]; // array containing all atoms in the simulation
double elapsed_time = 0;

bool evaporation_flag = true;
char coordinate_log_prefix[256] = "default_simulation_analysis.dat";

double run_time = 1.e8; //default (in seconds) // TODO: move to input file
double data_time_interval = 0.1;
double time_interval_end;

int rate_cnt; // number of rates in rate_arr list (filled indices)
int transition_cnt; // size of filled portion of transition list
double frequency_sum;

double overpotential = 0.0;

double nnE[6] = {
	DEFAULT_BOND_ENERGY_AA,
	DEFAULT_BOND_ENERGY_AB,
	DEFAULT_BOND_ENERGY_AC,
	DEFAULT_BOND_ENERGY_BB,
	DEFAULT_BOND_ENERGY_BC,
	DEFAULT_BOND_ENERGY_CC};

double nnnE[6] = {0., 0., 0., 0., 0., 0.};

bool solubility[3] = {false, false, false}; // whether ABC-type atoms can be dissolved/evaporated - all elements cannot dissolve by default

double temperature = DEFAULT_TEMPERATURE;

int dissolution = DISSOLUTION;

int final_config_neighbor_cnt;
int intial_config_neighbor_cnt;

// ENHANCE: malloc?
Rate rate_arr[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS]; // array containing all the unique rate constants and count of atoms that have that k and indices in transition_arr
// Rates default initalized to all zeros values?
Transition *transition_arr[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS]; // array of the possible atom transitions (atom index + jump offset index)
//Transition transition_arr[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
// contains all the same information in atom_arr[i].transition_indices
Trans_Prob transition_probability;


long int final_iteration = 1e9; // TODO: move to input file
int lastxt, lastyt, lastzt;

//bool simulation_is_going = false; //probably don't need this

double sum_of_rate_populations;
double current_probability;

// ENHANCE: pass struct with all simulation parameters as argument
unsigned long perform_simulation(void) //potentially FILE* as arguments
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

	int framenum = 0;

	elapsed_time = 0.0;
	//writes data time intervals and run time in original code

	//simulation_is_going = true;

	//print system and zone sizes to the file?

	/*if (num_sims > 0) //this probably will not happen
	{
		do_initialize_simulation(simulation_type);
		elapsed_time = 0.0;
	}*/ 

	// initialize simulation kinetics, and draw a picture

	//printf("transition time!\n");
	for (i=0;i<atom_cnt;++i)	
		refresh_transitions(i);			// resets all kinetic paramters

	organize(atom_arr, atom_cnt); //replacement for copy_xyz_to_coord but might not be necessary

	//printf("transitioned and organized\n");

	total_volume_dissolved = 0; //do i care

	calculate_internal_energy(atom_cnt);
	//printf("energy calculated\n");
	output_log_file(sim_log_file, framenum);
	write_xyz_file(coordinate_log_prefix, framenum);
	//printf("files written\n");
	++framenum;

	long int iter = 0;

	ot = 0.0; //needs to happen outside of loop
	//printf("about to start the loop\n");
			
	while (elapsed_time < run_time && iter < final_iteration) //adjusted this condition, included sanity check
	{
		if (iter % 100 == 0)
			printf("iteration %ld, time %lf\n", iter, elapsed_time);
		if (simulation_should_kill_itself)							// abort simulation (only happens if atoms overlap)
		{
			//find_average_curvature(); //no longer valid

			calculate_internal_energy(atom_cnt);
			output_log_file(sim_log_file, framenum);
			write_xyz_file(coordinate_log_prefix, framenum);

			simulation_should_kill_itself = false;
			organize(atom_arr, atom_cnt); //replace the copy with draw, maybe not needed

			//simulation_is_going = false; //not needed
			return 1; //return 1 b/c error?
		}


		for (j=0;j<atom_cnt;++j)	
			refresh_transitions(j);			// resets all kinetic paramters
		 //does this need to happen?

		// increment the elapsed time

		compute_transition_array();
		elapsed_time -= log(drandj(&rand_seed))/frequency_sum;

		if (elapsed_time >= time_interval_end) // [ ]: time interval for what???
		{
			organize(atom_arr, atom_cnt); //replaced but do i really need it
			printf("writing file %d: elapsed_time = %lf\n", framenum, elapsed_time);
			if (analysis_type == REGULAR_TIME_INTERVALS)
			{
				//record the elapsed time in a file here

				calculate_internal_energy(atom_cnt);
				output_log_file(sim_log_file, framenum);
				write_xyz_file(coordinate_log_prefix, framenum);

				while (time_interval_end <= elapsed_time)
					time_interval_end += data_time_interval;

				//write something with data_time_interval here?

				if (overpotential_ramp_rate != 0.0)
				{
					nt = elapsed_time;
					overpotential += (nt-ot)*overpotential_ramp_rate;
					ot = elapsed_time;
					for (j=0;j<atom_cnt;++j)	
						refresh_transitions(j);			// resets all kinetic paramters
	
				}
			}
			else if (analysis_type == LOG_TIME_INTERVALS)
			{
				calculate_internal_energy(atom_cnt);
				output_log_file(sim_log_file, framenum);
				write_xyz_file(coordinate_log_prefix, framenum);

				while (time_interval_end <= elapsed_time)
					time_interval_end *= logtime_multiplier;
			}
			++framenum;
		}

		if (elapsed_time >= run_time) // simulation has gone past time
			break; //get outta here before I make a new transition

		// pick the type of transition to occur
		transition_type_probability = drandj(&rand_seed);
		rate_skip = rate_cnt/2;
		j = rate_skip;

		//change structure to see if the move gets made or not

		// binary search to select transition
		moved_flag = false; //used to track diffusion/evaporation vs deposition
		while (moved_flag == false) {
			if ((transition_type_probability >= transition_probability.lbound[j])
						&& (transition_type_probability < transition_probability.ubound[j]))
			{
				k = transition_probability.listnum[j];

				// pick the lucky atom
				// [ ]: third random number?
				// picks a type of transition and then which atom that has that transition will it act on?
				which_one = rate_arr[k].offset + (int)(drandj(&rand_seed)*(double)rate_arr[k].count);

				// which_one gives the location of the transition_arr, which gives
				// the info about the specific atom
				atom_number = transition_arr[which_one]->atom_idx;
				jump_vector = transition_arr[which_one]->offset_idx;

				// if jump_vector == number_of_possible_neighbors then the atom is going to evaporate
				moved_flag = true;

				adatom_before = 0;

				if (jump_vector != number_of_possible_neighbors) //diffusion
				{
					//printf("the transition is diffusion of atom %d, jumping from %lf %lf %lf ", atom_number, atom[atom_number]->lattice[0], atom[atom_number]->lattice[1], atom[atom_number]->lattice[2]);
					// coordinates atom is jumping to
					lastxt = atom_arr[atom_number]->lattice[0] + jump_offset[jump_vector].dx;
					lastyt = atom_arr[atom_number]->lattice[1] + jump_offset[jump_vector].dy;
					lastzt = atom_arr[atom_number]->lattice[2] + jump_offset[jump_vector].dz;

					adjust_pbc(&lastxt, &lastyt, &lastzt);

					atype = atom_arr[atom_number]->type;

					// moves atom?
					remove_atom(atom_number);
					natn = add_atom(lastxt, lastyt, lastzt, atype, NORMAL);
					//printf("and jumping to %lf %lf %lf\n", lastxt, lastyt, lastzt);
				}
				else				// dissolution
				{
					if (solubility[atom_arr[atom_number]->type - 1] == true)	// atoms are dissolved based on input specs!
					{
						++total_volume_dissolved;
						remove_atom(atom_number);	// evaporate the atom
						//printf("the transition was dissolution of atom %d\n", atom_number);
					}
				}
			}
			else if (transition_type_probability < transition_probability.lbound[j])
			{
				//search to the left
				rate_skip = rate_skip/2;
				if (rate_skip == 0) rate_skip = 1;
				j -= rate_skip;
			}
			else if (transition_type_probability >= transition_probability.ubound[j])
			{
				//search to the right
				rate_skip = rate_skip/2;
				if (rate_skip == 0) rate_skip = 1;
				j += rate_skip;
				if (j == rate_cnt) //no more options!
					break;
			}
		}

		if (moved_flag == false) //only happens iff jump_vector == number_of_possible_neighbors
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

		++iter; //sanity check to avoid ending in an infinite cycle
	}
		if (iter == final_iteration)
			fprintf(sim_log_file, "reached final iteration and terminated\n");
		//write elapsed_time to mark finish


		//TODO: finish IO
		calculate_internal_energy(atom_cnt);
		output_log_file(sim_log_file, framenum);
		write_xyz_file(coordinate_log_prefix, framenum);	

	printf("Finished simulation\n"); //move this to the log file
		
	//simulation_is_going = false;
	return 0;
}

/******************************************************************************/
/******************************************************************************/

void compute_transition_array(void)
{
	int i;
	int x;

	int total_lists = 0;

	frequency_sum = 0.0;
	sum_of_rate_populations = 0.0;

	for (i=0;i<rate_cnt;++i)
	{
		if (rate_arr[i].count != 0)
		{
			rate_arr[i].frequency = rate_arr[i].k*(double)rate_arr[i].count;
			frequency_sum += rate_arr[i].frequency;
			sum_of_rate_populations += rate_arr[i].count;

			transition_probability.listnum[total_lists] = i;
			++total_lists;
		}
	}

	/*if (deposition_type == DEPOSITION_TYPE_RAINFALL)
		frequency_sum += (double)ssx*(double)ssy*(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c);
	else if (deposition_type == DEPOSITION_TYPE_RANDOM_WALKER)
		frequency_sum += (4*PI*ssr*ssr)*(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c);*/ //deposition is vestigial

	// now compute bounds for jump probabilities

	current_probability = 0.0;

	for (i=0;i<total_lists;++i)
	{
		transition_probability.lbound[i] = current_probability;
		x = transition_probability.listnum[i];
		current_probability += rate_arr[x].frequency/frequency_sum;
		transition_probability.ubound[i] = current_probability;
	}

	return;
}

/******************************************************************************/
/******************************************************************************/
// updates [atom_arr[atom_idx], transition_arr[i], rate_arr[i].offset]
int refresh_transitions(int atom_idx) // atom_idx = index on atom list
{
	int i, j;
	int rate_idx; // position of rate rate in rate list rate_arr[]
	double rate; // rate constant from bond-breaking model
	int next_x, next_y, next_z;

	int atom_rates_cnt;	// this is returned as the number of transitions (jump_offsets) this atom can undergo, excluding evaporation
	// 
	int start_config[MAXIMUM_NUMBER_OF_NEIGHBORS]; // -1 if empty, type if filled
	int end_config[MAXIMUM_NUMBER_OF_NEIGHBORS];

	//bool ok_to_evaporate; //doesn't get used // XXX: commented code

	//printf("removing\n"); // XXX: commented print
	// first, remove all mention of this atom from transition list
	for (i=0; i<number_of_possible_neighbors + dissolution; ++i) // extra 1 for evaporation
	{
		if (atom_arr[atom_idx]->transition_indices[i] != -1) // something can happen in the "i" direction
			take_off_transition_list(atom_idx, i);
	}

	//printf("cycling\n"); // XXX: commented print
	// cycle through neighbor coordinates, and check if there is an atom there
	for (i=0;i<number_of_possible_neighbors;++i)
	{
   		next_x = atom_arr[atom_idx]->lattice[0] + jump_offset[i].dx;
		next_y = atom_arr[atom_idx]->lattice[1] + jump_offset[i].dy;
        next_z = atom_arr[atom_idx]->lattice[2] + jump_offset[i].dz;

		adjust_pbc(&next_x, &next_y, &next_z);

        j = atom_at(next_x, next_y, next_z);

		if (j >= 0) atom_arr[atom_idx]->occupied_neighbor_sites[i] = j;
		// if no atom atom_idx that position but occ_neighbor array says there is, fix it
		if ((j == -1)&&(atom_arr[atom_idx]->occupied_neighbor_sites[i] >= 0))
			atom_arr[atom_idx]->occupied_neighbor_sites[i] = -1;
	}

	// cycle through the neighbor sites.  if there's an empty one, calculate the transition rate to it

	//printf("getting config\n"); // XXX: commented prints
	atom_rates_cnt = 0;
	//ok_to_evaporate = true; //doesn't get used

	intial_config_neighbor_cnt = get_initial_configuration2(atom_idx, 0, start_config);		// k is number of near neighbors

	if (intial_config_neighbor_cnt == number_of_possible_neighbors) // skip calculating a rate of a fully coordinated atom
		return atom_rates_cnt;			

	//printf("calculating surf diffusion rate\n");
	for (i=0; i<number_of_possible_neighbors; ++i)
	{
		if (atom_arr[atom_idx]->occupied_neighbor_sites[i] == -1)
		{ // if unoccupied, consider transition
			//printf("final config\n");
			final_config_neighbor_cnt = get_final_configuration2(atom_idx, i, end_config);
			//printf("surf diffusion\n");
			calculate_surf_diffusion_rate(		start_config, // ENHANCE: make this look prettier
												end_config,
												number_of_possible_neighbors,
												atom_arr[atom_idx]->type,
												nnE,
												temperature,
												overpotential,
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
								number_of_possible_neighbors,
								atom_arr[atom_idx]->type,
								nnE,
								temperature,
								overpotential,
								&rate);								

	//replace with the lines below because it's more obvious
	/*if ((rate_idx = is_on_transition_list(rate)) != -1)
	{
		add_to_transition_list(rate_idx, atom_idx, number_of_possible_neighbors);
	}
	else  // the transition rate to that spot turned out to be a new one.
	{		
		rate_idx = create_new_transition(rate);
		add_to_transition_list(rate_idx, atom_idx, number_of_possible_neighbors);
	}*/

	//printf("end stuff\n");
	rate_idx = is_on_transition_list(rate);

	if (rate_idx == -1)
		rate_idx = create_new_transition(rate); //the transition rate to the spot is a new one!

	add_to_transition_list(rate_idx, atom_idx, number_of_possible_neighbors); // evaporation is considered to be last in jump_offset (not really in array but uses that index number)
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

	for (i=0;i<rate_cnt;++i)
		if (rate == rate_arr[i].k) return i;

	return -1;
}

/******************************************************************************/
/******************************************************************************/
// create new Rate struct in rate array
// updates rate_array[rate_cnt], rate_cnt
int create_new_transition(double rate)
{
	rate_arr[rate_cnt].k = rate;
	rate_arr[rate_cnt].offset = transition_cnt;
	rate_arr[rate_cnt].count = 0;

	++rate_cnt;

	return (rate_cnt-1);
}

/******************************************************************************/
/******************************************************************************/

// add to rate_arr[rate_idx] the atom atom_idx going in direction offset_idx
// updates transition_arr, rate_arr[rate_idx].count, atom_arr[atom_idx]->transition_indices[offset_idx]
void add_to_transition_list(int rate_idx, int atom_idx, int offset_idx) // rate_arr index, atom_arr index, jump_offset index
{
	int i; // loop variable
	int n;
	int initial_transition_index, final_transition_index; // initial and final transition_arr index

	// make room for the new arrival

	transition_arr[transition_cnt] = (Transition *)malloc(sizeof(Transition));	// adds entry to the end of the list

	// what is this
	//		final_transition_index = rate_arr[rate_cnt-1].offset + rate_arr[rate_cnt-1].number;	
	//		transition_arr[final_transition_index] = (Transition *)malloc(sizeof(Transition));

	for (i = rate_cnt-1;i>rate_idx;--i)
		{ // [ ]: what does this do? is this the same as in remove_transition?
			initial_transition_index = rate_arr[i].offset;
			final_transition_index = initial_transition_index + rate_arr[i].count;

			transition_arr[final_transition_index]->atom_idx = transition_arr[initial_transition_index]->atom_idx;
			transition_arr[final_transition_index]->offset_idx = transition_arr[initial_transition_index]->offset_idx;

			if (initial_transition_index != final_transition_index)
				atom_arr[transition_arr[final_transition_index]->atom_idx]->transition_indices[transition_arr[final_transition_index]->offset_idx] = final_transition_index;

			++rate_arr[i].offset;
		}

	// add new arrival

	n = rate_arr[rate_idx].offset + rate_arr[rate_idx].count;

	++rate_arr[rate_idx].count;

	transition_arr[n]->atom_idx = atom_idx;
	transition_arr[n]->offset_idx = offset_idx;

	atom_arr[atom_idx]->transition_indices[offset_idx] = n;

	++transition_cnt;

	return;
}

/******************************************************************************/
/******************************************************************************/
// updates atom_arr[atom_idx], transition_arr[i], rate_arr[i].offset
void take_off_transition_list(int atom_idx, int offset_idx)	// removes atom jumping in the jump_offset[offset_idx] direction
{
	int i;
	int rate_idx, transition_start_idx, transition_end_idx;

	// find out what Rate in rate_arr this is

	transition_start_idx = atom_arr[atom_idx]->transition_indices[offset_idx];

	for (i=0;i<rate_cnt;++i)
		if (transition_start_idx < (rate_arr[i].offset + rate_arr[i].count))
		{
			rate_idx = i;
			break;
		}

	// remind atom it can no longer jump

	atom_arr[atom_idx]->transition_indices[offset_idx] = -1;
	// TODO: do these later, after it's been removed and transition_arr has been rearranged
	--transition_cnt;

	// rate_idx points to the current rate list it's on.  decrement the number of atoms in that list
	// and clean up.  If the list is empty, remove it.

	--rate_arr[rate_idx].count;
	// [ ]: wtf is going on here 
	if (rate_arr[rate_idx].count == 0) // if list is empty
	{
		for (i = rate_idx + 1; i < rate_cnt; ++i)
		{
			--rate_arr[i].offset; // move rate_arr offsets of larger indicies down one
			// make transition at new start index (which is of different Rate than old start index) have same atom and offset as new end index (which is of same Rate as old end index)
			transition_start_idx = rate_arr[i].offset;
			transition_end_idx = rate_arr[i].offset+rate_arr[i].count;		// count is always at least 1

			transition_arr[transition_start_idx]->atom_idx = transition_arr[transition_end_idx]->atom_idx;
			transition_arr[transition_start_idx]->offset_idx = transition_arr[transition_end_idx]->offset_idx;
			
			atom_arr[transition_arr[transition_start_idx]->atom_idx]->transition_indices[transition_arr[transition_start_idx]->offset_idx] = transition_start_idx;
		}

		free(transition_arr[transition_cnt]);			// free up the very last member of the last transition_arr

		for (i=rate_idx+1;i<rate_cnt;++i)
		{
			rate_arr[i-1].offset = rate_arr[i].offset;
			rate_arr[i-1].count = rate_arr[i].count;
			rate_arr[i-1].k = rate_arr[i].k;
			rate_arr[i-1].frequency = rate_arr[i].frequency;
		}

		--rate_cnt;

		return;
	}

	transition_end_idx = rate_arr[rate_idx].offset + rate_arr[rate_idx].count; // points to what was last element

	// swap transition_end_idx into the position atom_idx:offset_idx occupied
	// ENHANCE: this is the same shit that happens when count==0
	transition_arr[transition_start_idx]->atom_idx = transition_arr[transition_end_idx]->atom_idx;
	transition_arr[transition_start_idx]->offset_idx = transition_arr[transition_end_idx]->offset_idx;

	if (transition_start_idx != transition_end_idx) 
		atom_arr[transition_arr[transition_start_idx]->atom_idx]->transition_indices[transition_arr[transition_start_idx]->offset_idx] = transition_start_idx;

	// shift all other transition lists

	for (i=rate_idx+1;i < rate_cnt;++i)
		{ // ENHANCE: again, looks like the same shit that happens when count==0
			--rate_arr[i].offset;

			transition_start_idx = rate_arr[i].offset;
			transition_end_idx = rate_arr[i].offset+rate_arr[i].count;

			transition_arr[transition_start_idx]->atom_idx = transition_arr[transition_end_idx]->atom_idx;
			transition_arr[transition_start_idx]->offset_idx = transition_arr[transition_end_idx]->offset_idx;

			atom_arr[transition_arr[transition_start_idx]->atom_idx]->transition_indices[transition_arr[transition_start_idx]->offset_idx] = transition_start_idx;
		}

	free(transition_arr[transition_cnt]);			// free up the very last member of the last transition_arr

	return;
}


/******************************************************************************/
/******************************************************************************/

void check_system(void)
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

	for (j=0;j<atom_cnt;++j)
		for (i=0;i<number_of_possible_neighbors + 1;++i)				// extra 1 for evaporation
			{
				if (atom_arr[j]->transition_indices[i] != -1)		// something can happen in the "i" direction
					take_off_transition_list(j, i);
			}

	// for each atom, cycle through neighbor coordinates, and and reconcile occupancy

	for (j=0;j<atom_cnt;++j)
		for (i=0;i<number_of_possible_neighbors;++i)
		{
   			next_x = atom_arr[j]->lattice[0] + jump_offset[i].dx;
			next_y = atom_arr[j]->lattice[1] + jump_offset[i].dy;
	        next_z = atom_arr[j]->lattice[2] + jump_offset[i].dz;

			adjust_pbc(&next_x, &next_y, &next_z);

			k = atom_at(next_x, next_y, next_z);

			if (k >= 0)
			{
				// an atom has been found atom_idx this neighbor site.

				atom_arr[j]->occupied_neighbor_sites[i] = k;
				atom_arr[k]->occupied_neighbor_sites[opposite_offset[i]] = j;
			}
		}

	// now we'll reconcile buried atoms.
	// does this process need to happen if nothing is buried???
	do
	{
		errors = 0;

		for (j=0;j<atom_cnt;++j)
		for (i=0;i<number_of_possible_neighbors;++i)
		{
			if (atom_arr[j]->occupied_neighbor_sites[i] == -2)
			{
				// find coordinate of buried atom

				next_x = atom_arr[j]->lattice[0] + jump_offset[i].dx;
				next_y = atom_arr[j]->lattice[1] + jump_offset[i].dy;
		        next_z = atom_arr[j]->lattice[2] + jump_offset[i].dz;

				adjust_pbc(&next_x, &next_y, &next_z);

				for (k=0;k<number_of_possible_neighbors;++k)
				{
					nnx = next_x + jump_offset[k].dx;
					nny = next_y + jump_offset[k].dy;
					nnz = next_z + jump_offset[k].dz;

					adjust_pbc(&nnx, &nny, &nnz);

					m = atom_at(nnx, nny, nnz);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]];
						if (n != -2)
						{
							if (n == -1) atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]] = -2;
							if (n == -3) atom_arr[j]->occupied_neighbor_sites[i] = -3;		// random trumps buried

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

		for (j=0;j<atom_cnt;++j)
		for (i=0;i<number_of_possible_neighbors;++i)
		{
			if (atom_arr[j]->occupied_neighbor_sites[i] == -3)
			{
				// find coordinate of buried atom

				next_x = atom_arr[j]->lattice[0] + jump_offset[i].dx;
				next_y = atom_arr[j]->lattice[1] + jump_offset[i].dy;
		        next_z = atom_arr[j]->lattice[2] + jump_offset[i].dz;

				adjust_pbc(&next_x, &next_y, &next_z);

				for (k=0;k<number_of_possible_neighbors;++k)
				{
					nnx = next_x + jump_offset[k].dx;
					nny = next_y + jump_offset[k].dy;
					nnz = next_z + jump_offset[k].dz;

					adjust_pbc(&nnx, &nny, &nnz);

					m = atom_at(nnx, nny, nnz);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]];
						if (n != -3)
						{
							atom_arr[m]->occupied_neighbor_sites[opposite_offset[k]] = -3;
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

		for (i=0;i<number_of_possible_neighbors;++i)
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

		if (k == number_of_possible_neighbors)
		{
			// bury atom j

			for (i=0;i<number_of_possible_neighbors;++i)
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

	for (j=0;j<atom_cnt;++j)
		refresh_transitions(j);

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
	for (i=0;i<number_of_possible_neighbors;++i)
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
	//*rate = 1e13*exp(-energy/(kBoltz*temperature)) //+ 1e-4*exp(-(energy-overpotential)/(kBoltz*temperature));
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

	/*if ((ncsk == 1)||(evaporation_flag == false)) //I don't think this happens
	{
		energy = 1000.;
		*rate = 1e4*exp(-(energy)/(kBoltz*temperature));
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
			if (solubility[0]) {
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
			if (solubility[1]) {
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
			if (solubility[2]) {
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