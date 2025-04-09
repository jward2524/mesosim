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
int nat = 0;
Atom temp_atom;
bool simulation_should_kill_itself;
Atom* atom[];
double elapsed_time = 0;

bool evaporation_flag = true;
char coordinate_log_prefix[256] = "default_simulation_analysis.dat";

double run_time = 1.e8; //default (in seconds)
double data_time_interval = 0.1;
double time_interval_end;

int number_rates;
int total_current_transitions;
double sum_of_frequencies;

double overpotential = 0.0;

double nnE[6] = {
	DEFAULT_BOND_ENERGY_AA,
	DEFAULT_BOND_ENERGY_AB,
	DEFAULT_BOND_ENERGY_AC,
	DEFAULT_BOND_ENERGY_BB,
	DEFAULT_BOND_ENERGY_BC,
	DEFAULT_BOND_ENERGY_CC};

double nnnE[6] = {0., 0., 0., 0., 0., 0.};

bool solubility[3] = {false, false, false}; //all elements cannot dissolve by default

double temperature = DEFAULT_TEMPERATURE;

int dissolution = DISSOLUTION;

int number_final_configuration_neighbors;
int number_intial_configuration_neighbors;

// ENHANCE: malloc?
Rate rate[MAXIMUM_NUMBER_OF_ACTIVATION_BARRIERS];

Transition_List *transition_list[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];
//Transition_List transition_list[MAXIMUM_NUMBER_OF_CONCURRENT_TRANSITIONS];

Trans_Prob transition_probability;


long int final_iteration = 1e9;
double lastxt, lastyt, lastzt;

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
	for (i=0;i<nat;++i)	
		refresh_transitions(i);			// resets all kinetic paramters

	organize(atom, nat); //replacement for copy_xyz_to_coord but might not be necessary

	//printf("transitioned and organized\n");

	total_volume_dissolved = 0; //do i care

	calculate_internal_energy(nat);
	//printf("energy calculated\n");
	output_log_file(framenum);
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

			calculate_internal_energy(nat);
			output_log_file(framenum);
			write_xyz_file(coordinate_log_prefix, framenum);

			simulation_should_kill_itself = false;
			organize(atom, nat); //replace the copy with draw, maybe not needed

			//simulation_is_going = false; //not needed
			return 1; //return 1 b/c error?
		}


		for (j=0;j<nat;++j)	
			refresh_transitions(j);			// resets all kinetic paramters
		 //does this need to happen?

		// increment the elapsed time

		compute_transition_array();
		elapsed_time -= log(drandj(&seed))/sum_of_frequencies;

		if (elapsed_time >= time_interval_end) // [ ]: time interval for what???
		{
			organize(atom, nat); //replaced but do i really need it
			printf("writing file %d: elapsed_time = %lf\n", framenum, elapsed_time);
			if (analysis_type == REGULAR_TIME_INTERVALS)
			{
				//record the elapsed time in a file here

				calculate_internal_energy(nat);
				output_log_file(framenum);
				write_xyz_file(coordinate_log_prefix, framenum);

				while (time_interval_end <= elapsed_time)
					time_interval_end += data_time_interval;

				//write something with data_time_interval here?

				if (overpotential_ramp_rate != 0.0)
				{
					nt = elapsed_time;
					overpotential += (nt-ot)*overpotential_ramp_rate;
					ot = elapsed_time;
					for (j=0;j<nat;++j)	
						refresh_transitions(j);			// resets all kinetic paramters
	
				}
			}
			else if (analysis_type == LOG_TIME_INTERVALS)
			{
				calculate_internal_energy(nat);
				output_log_file(framenum);
				write_xyz_file(coordinate_log_prefix, framenum);

				while (time_interval_end <= elapsed_time)
					time_interval_end *= logtime_multiplier;
			}
			++framenum;
		}

		if (elapsed_time >= run_time) // simulation has gone past time
			break; //get outta here before I make a new transition

		// pick the type of transition to occur
		transition_type_probability = drandj(&seed);
		rate_skip = number_rates/2;
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
				which_one = rate[k].offset + (int)(drandj(&seed)*(double)rate[k].number);

				// which_one gives the location of the transition_list, which gives
				// the info about the specific atom
				atom_number = transition_list[which_one]->number_in_list;
				jump_vector = transition_list[which_one]->offset_vector;

				// if jump_vector == number_of_possible_neighbors then the atom is going to evaporate
				moved_flag = true;

				adatom_before = 0;

				if (jump_vector != number_of_possible_neighbors) //diffusion
				{
					//printf("the transition is diffusion of atom %d, jumping from %lf %lf %lf ", atom_number, atom[atom_number]->lattice[0], atom[atom_number]->lattice[1], atom[atom_number]->lattice[2]);
					// coordinates atom is jumping to
					lastxt = atom[atom_number]->lattice[0] + jump_offset[jump_vector].dx;
					lastyt = atom[atom_number]->lattice[1] + jump_offset[jump_vector].dy;
					lastzt = atom[atom_number]->lattice[2] + jump_offset[jump_vector].dz;

					adjust_pbc(&lastxt, &lastyt, &lastzt);

					atype = atom[atom_number]->type;

					// moves atom?
					remove_atom(atom_number);
					natn = add_atom(lastxt, lastyt, lastzt, atype, NORMAL);
					//printf("and jumping to %lf %lf %lf\n", lastxt, lastyt, lastzt);
				}
				else				// dissolution
				{
					if (solubility[atom[atom_number]->type - 1] == true)	// atoms are dissolved based on input specs!
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
				if (j == number_rates) //no more options!
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

				vap = drandj(&seed);
				if (vap < deposition_rate_of_a/(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c))
					atype = 1; 
				else if (vap < (deposition_rate_of_a + deposition_rate_of_b)/(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c))
					atype = 2;
				else atype = 3;

				natn = add_atom(dep_x, dep_y, dep_z, atype, NORMAL);
			}
			else if (deposition_type == DEPOSITION_TYPE_RANDOM_WALKER)
			{
				vap = drandj(&seed);
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
					ta1 = 2.*PI*drandj(&seed);
					ta2 = 2.*PI*drandj(&seed);

					rrp[0] = ssr*cos(ta2)*cos(ta1);
					rrp[1] = ssr*cos(ta2)*sin(ta1);
					rrp[2] = ssr*sin(ta2);

					// now invert xr, yr, zr into lattice vectors;

					vecmul(rrp, ilatmat, rw);

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
		calculate_internal_energy(nat);
		output_log_file(framenum);
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

	sum_of_frequencies = 0.0;
	sum_of_rate_populations = 0.0;

	for (i=0;i<number_rates;++i)
	{
		if (rate[i].number != 0)
		{
			rate[i].frequency = rate[i].k*(double)rate[i].number;
			sum_of_frequencies += rate[i].frequency;
			sum_of_rate_populations += rate[i].number;

			transition_probability.listnum[total_lists] = i;
			++total_lists;
		}
	}

	/*if (deposition_type == DEPOSITION_TYPE_RAINFALL)
		sum_of_frequencies += (double)ssx*(double)ssy*(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c);
	else if (deposition_type == DEPOSITION_TYPE_RANDOM_WALKER)
		sum_of_frequencies += (4*PI*ssr*ssr)*(deposition_rate_of_a + deposition_rate_of_b + deposition_rate_of_c);*/ //deposition is vestigial

	// now compute bounds for jump probabilities

	current_probability = 0.0;

	for (i=0;i<total_lists;++i)
	{
		transition_probability.lbound[i] = current_probability;
		x = transition_probability.listnum[i];
		current_probability += rate[x].frequency/sum_of_frequencies;
		transition_probability.ubound[i] = current_probability;
	}

	return;
}

/******************************************************************************/
/******************************************************************************/

int refresh_transitions(int at)
{
	int i, j;
	int tl;
	double r;
	double nx, ny, nz;

	int nr;				// this is returned as the number of transitions this atom can undergo
		
	int initial_config[MAXIMUM_NUMBER_OF_NEIGHBORS]; // -1 if empty, type if filled
	int final_config[MAXIMUM_NUMBER_OF_NEIGHBORS];

	//bool ok_to_evaporate; //doesn't get used

	//printf("removing\n");
	// first, remove all mention of this atom from transition list
	for (i=0;i<number_of_possible_neighbors + dissolution;++i)		// extra 1 for evaporation
	{
		if (atom[at]->position_on_transition_list[i] != -1)		// something can happen in the "i" direction
			take_off_transition_list(at, i);
	}

	//printf("cycling\n");
	// cycle through neighbor coordinates, and check if there is an atom there
	for (i=0;i<number_of_possible_neighbors;++i)
	{
   		nx = atom[at]->lattice[0] + jump_offset[i].dx;
		ny = atom[at]->lattice[1] + jump_offset[i].dy;
        nz = atom[at]->lattice[2] + jump_offset[i].dz;

		adjust_pbc(&nx, &ny, &nz);

        j = atom_at(nx, ny, nz);

		if (j >= 0) atom[at]->occupied_neighbor_sites[i] = j;

		if ((j == -1)&&(atom[at]->occupied_neighbor_sites[i] >= 0))
			atom[at]->occupied_neighbor_sites[i] = -1;
	}

	// cycle through the neighbor sites.  if there's an empty one, calculate the transition rate to it

	//printf("getting config\n");
	nr = 0;
	//ok_to_evaporate = true; //doesn't get used

	number_intial_configuration_neighbors = get_initial_configuration2(at, 0, initial_config);		// k is number of near neighbors

	if (number_intial_configuration_neighbors == number_of_possible_neighbors) return nr;			// no sense in calculating a rate of a fully coordinated atom

	//printf("calculating surf diffusion rate\n");
	for (i=0;i<number_of_possible_neighbors;++i)
	{
		if (atom[at]->occupied_neighbor_sites[i] == -1)
		{
			//printf("final config\n");
			number_final_configuration_neighbors = get_final_configuration2(at, i, final_config);
			//printf("surf diffusion\n");
			calculate_surf_diffusion_rate(		initial_config,
												final_config,
												number_of_possible_neighbors,
												atom[at]->type,
												nnE,
												temperature,
												overpotential,
												&r);
					
			++nr;
			//printf("adding step\n");
			if ((tl = is_on_transition_list(r)) != -1)
			{
				add_to_transition_list(tl, at, i);
			}
			else
			{
				// the transition rate to that spot turned out to be a new one.

				tl = create_new_transition(r);
				add_to_transition_list(tl, at, i);
			}
		}
	}	

	//printf("calcualte evap\n");
	calculate_evaporation_rate(	initial_config,
								number_of_possible_neighbors,
								atom[at]->type,
								nnE,
								temperature,
								overpotential,
								&r);								

	//replace with the lines below because it's more obvious
	/*if ((tl = is_on_transition_list(r)) != -1)
	{
		add_to_transition_list(tl, at, number_of_possible_neighbors);
	}
	else  // the transition rate to that spot turned out to be a new one.
	{		
		tl = create_new_transition(r);
		add_to_transition_list(tl, at, number_of_possible_neighbors);
	}*/

	//printf("end stuff\n");
	tl = is_on_transition_list(r);

	if (tl == -1)
		tl = create_new_transition(r); //the transition rate to the spot is a new one!

	add_to_transition_list(tl, at, number_of_possible_neighbors);
	//printf("gonna return %d\n", nr);
	return nr;						// gives number of current transitions for that atom
}

/******************************************************************************/
/******************************************************************************/

//is there a better way of doing this?

int is_on_transition_list(double r) 
{
	int i;

	for (i=0;i<number_rates;++i)
		if (r == rate[i].k) return i;

	return -1;
}

/******************************************************************************/
/******************************************************************************/

int create_new_transition(double r)
{
	rate[number_rates].k = r;
	rate[number_rates].offset = total_current_transitions;
	rate[number_rates].number = 0;

	++number_rates;

	return (number_rates-1);
}

/******************************************************************************/
/******************************************************************************/

// add to list[tl] the atom at going in direction vc
// [ ]: help understanding this - what is offset, number, ia, fa
void add_to_transition_list(int tl, int at, int vc) // transition list index, atom index, diffusion vector index
{
	int i;
	int n;
	int ia, fa;

	// make room for the new arrival

	transition_list[total_current_transitions] = (Transition_List *)malloc(sizeof(Transition_List));	// adds entry to the end of the list

	// what is this
	//		fa = rate[number_rates-1].offset + rate[number_rates-1].number;	
	//		transition_list[fa] = (Transition_List *)malloc(sizeof(Transition_List));

	for (i = number_rates-1;i>tl;--i)
		{
			ia = rate[i].offset;
			fa = ia+rate[i].number;

			transition_list[fa]->number_in_list = transition_list[ia]->number_in_list;
			transition_list[fa]->offset_vector = transition_list[ia]->offset_vector;

			if (ia != fa)
				atom[transition_list[fa]->number_in_list]->position_on_transition_list[transition_list[fa]->offset_vector] = fa;

			++rate[i].offset;
		}

	// add new arrival

	n = rate[tl].offset + rate[tl].number;

	++rate[tl].number;

	transition_list[n]->number_in_list = at;
	transition_list[n]->offset_vector = vc;

	atom[at]->position_on_transition_list[vc] = n;

	++total_current_transitions;

	return;
}

/******************************************************************************/
/******************************************************************************/

void take_off_transition_list(int at, int vc)			// removes atom jumping in the vc direction
	{
		int i;
		int tl, tl_pos1, tl_pos2;

		// find out what transition list this is

		tl_pos1 = atom[at]->position_on_transition_list[vc];

		for (i=0;i<number_rates;++i)
			if (tl_pos1 < (rate[i].offset + rate[i].number))
				{
					tl = i;
					break;
				}

		// remind atom it can no longer jump

		atom[at]->position_on_transition_list[vc] = -1;

		--total_current_transitions;

		// tl points to the current list it's on.  decrement the number of atoms in that list
		// and clean up.  If the list is empty, remove it.

		--rate[tl].number;

		if (rate[tl].number == 0)
			{
				for (i=tl+1;i < number_rates;++i)
					{
						--rate[i].offset;

						tl_pos1 = rate[i].offset;
						tl_pos2 = rate[i].offset+rate[i].number;		// number is always at least one.

						transition_list[tl_pos1]->number_in_list = transition_list[tl_pos2]->number_in_list;
						transition_list[tl_pos1]->offset_vector = transition_list[tl_pos2]->offset_vector;

						atom[transition_list[tl_pos1]->number_in_list]->position_on_transition_list[transition_list[tl_pos1]->offset_vector] = tl_pos1;
					}

				free(transition_list[total_current_transitions]);			// free up the very last member of the last transition_list

				for (i=tl+1;i<number_rates;++i)
					{
						rate[i-1].offset = rate[i].offset;
						rate[i-1].number = rate[i].number;
						rate[i-1].k = rate[i].k;
						rate[i-1].frequency = rate[i].frequency;
					}

				--number_rates;

				return;
			}

		tl_pos2 = rate[tl].offset + rate[tl].number;		// points to what was last element

		// swap t1_pos2 into the position our atom occupied.

		transition_list[tl_pos1]->number_in_list =	transition_list[tl_pos2]->number_in_list;
		transition_list[tl_pos1]->offset_vector = transition_list[tl_pos2]->offset_vector;

		if (tl_pos1 != tl_pos2)
			atom[transition_list[tl_pos1]->number_in_list]->position_on_transition_list[transition_list[tl_pos1]->offset_vector] = tl_pos1;

		// shift all other transition lists

		for (i=tl+1;i < number_rates;++i)
			{
				--rate[i].offset;

				tl_pos1 = rate[i].offset;
				tl_pos2 = rate[i].offset+rate[i].number;

				transition_list[tl_pos1]->number_in_list = transition_list[tl_pos2]->number_in_list;
				transition_list[tl_pos1]->offset_vector = transition_list[tl_pos2]->offset_vector;

				atom[transition_list[tl_pos1]->number_in_list]->position_on_transition_list[transition_list[tl_pos1]->offset_vector] = tl_pos1;
			}

		free(transition_list[total_current_transitions]);			// free up the very last member of the last transition_list

		return;
	}


/******************************************************************************/
/******************************************************************************/

void check_system(void)
{
	int i, j, k, m, n, mm;
	int errors;

	double nx, ny, nz;

	double nnx, nny, nnz;

	int xzone, yzone, zzone;

	double xx1[3], xx2[3], xx3[3];
	double yy1[3], yy2[3], yy3[3];

	// does a careful check to make sure that a system is ready to be simulated.
	// assumptions:  (1) all real atoms are in the places they think they are
	// (2) all buried atoms are actually buried.

	// first, remove all atoms from the transition list.  We'll add them after we check neighbors

	for (j=0;j<nat;++j)
		for (i=0;i<number_of_possible_neighbors + 1;++i)				// extra 1 for evaporation
			{
				if (atom[j]->position_on_transition_list[i] != -1)		// something can happen in the "i" direction
					take_off_transition_list(j, i);
			}

	// for each atom, cycle through neighbor coordinates, and and reconcile occupancy

	for (j=0;j<nat;++j)
		for (i=0;i<number_of_possible_neighbors;++i)
		{
   			nx = atom[j]->lattice[0] + jump_offset[i].dx;
			ny = atom[j]->lattice[1] + jump_offset[i].dy;
	        nz = atom[j]->lattice[2] + jump_offset[i].dz;

			adjust_pbc(&nx, &ny, &nz);

			k = atom_at(nx, ny, nz);

			if (k >= 0)
			{
				// an atom has been found at this neighbor site.

				atom[j]->occupied_neighbor_sites[i] = k;
				atom[k]->occupied_neighbor_sites[opposite_offset[i]] = j;
			}
		}

	// now we'll reconcile buried atoms.
	// does this process need to happen if nothing is buried???
	do
	{
		errors = 0;

		for (j=0;j<nat;++j)
		for (i=0;i<number_of_possible_neighbors;++i)
		{
			if (atom[j]->occupied_neighbor_sites[i] == -2)
			{
				// find coordinate of buried atom

				nx = atom[j]->lattice[0] + jump_offset[i].dx;
				ny = atom[j]->lattice[1] + jump_offset[i].dy;
		        nz = atom[j]->lattice[2] + jump_offset[i].dz;

				adjust_pbc(&nx, &ny, &nz);

				for (k=0;k<number_of_possible_neighbors;++k)
				{
					nnx = nx + jump_offset[k].dx;
					nny = ny + jump_offset[k].dy;
					nnz = nz + jump_offset[k].dz;

					adjust_pbc(&nnx, &nny, &nnz);

					m = atom_at(nnx, nny, nnz);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = atom[m]->occupied_neighbor_sites[opposite_offset[k]];
						if (n != -2)
						{
							if (n == -1) atom[m]->occupied_neighbor_sites[opposite_offset[k]] = -2;
							if (n == -3) atom[j]->occupied_neighbor_sites[i] = -3;		// random trumps buried

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

		for (j=0;j<nat;++j)
		for (i=0;i<number_of_possible_neighbors;++i)
		{
			if (atom[j]->occupied_neighbor_sites[i] == -3)
			{
				// find coordinate of buried atom

				nx = atom[j]->lattice[0] + jump_offset[i].dx;
				ny = atom[j]->lattice[1] + jump_offset[i].dy;
		        nz = atom[j]->lattice[2] + jump_offset[i].dz;

				adjust_pbc(&nx, &ny, &nz);

				for (k=0;k<number_of_possible_neighbors;++k)
				{
					nnx = nx + jump_offset[k].dx;
					nny = ny + jump_offset[k].dy;
					nnz = nz + jump_offset[k].dz;

					adjust_pbc(&nnx, &nny, &nnz);

					m = atom_at(nnx, nny, nnz);

					if ((m >= 0)&&(m!= j))
					{
						// another atom (m) is connected to this atom.  If it sees this position as 
						// a buried atom, great.  Otherwise, reconcile

						n = atom[m]->occupied_neighbor_sites[opposite_offset[k]];
						if (n != -3)
						{
							atom[m]->occupied_neighbor_sites[opposite_offset[k]] = -3;
							++errors;
						}
					
					}
				}
			}
		}
	}
	while (errors != 0);

	// now let's bury any atoms that should be buried - DON'T WANT THIS NOW!
	/*for (j=0;j<nat;++j)
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

			if (j != (nat-1))
					move_atom((nat-1), j);

			free(atom[nat-1]);
			--nat;
			--j;
		}
	}*/

	// now we can recalculate diffusion rates

	for (j=0;j<nat;++j)
		refresh_transitions(j);

	return;
}


/******************************************************************************/
/******************************************************************************/

int calculate_surf_diffusion_rate(	int initial_configuration[],			// initial configuration array
									int final_configuration[],				// final configuration array
									int number_of_neighbors,				// number of 1st near neighbors in current xtal structure
									int atom_type,							// type of atom in consideration for transition
									double nnE[6],							// bond energy (type 2)-(type 2)
									double temperature,						// system temperature
									double overpotential,					// system overpotential
									double *rate)							// return value
{
	int i;
	double energy = 0.0;

	double nsi[3], nsf[3];
	double ti;
	double tf;

	static double baf[12] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};			// optional anistropy factor

	tf = 0.0;
	ti = 0.0;

	for (i=0;i<3;++i)
	{
		nsi[i] = 0.0;
		nsf[i] = 0.0;
	}

	int nneA_index[3] = {0, 1, 2}; //indices of A-A, A-B, A-C bonds
	int nneB_index[3] = {1, 3, 4}; //indices of B-A, B-B, B-C bonds
	int nneC_index[3] = {2, 4, 5}; //indices of C-A, C-B, C-C bonds
	int neighbor_type; //type of atom for nearest neighbor

	//printf("before teh switch: atom type = %d\n", atom_type);

	/*printf("final configuration: ");
	for (i=0;i<number_of_possible_neighbors;++i)
		printf("%d ", final_configuration[i]);
	printf("\n");*/

	// ENHANCE: lol these are basically identical, do it better (fxn)
	switch(atom_type)
	{
		case 1:
			for (i=0;i<number_of_neighbors;++i)
			{
				neighbor_type = initial_configuration[i];
				//printf("init neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++nsi[neighbor_type - 1]; //number of type A/B/C neighboring atoms in init configuration
					energy += nnE[nneA_index[neighbor_type - 1]]*baf[i]; //A-A/B/C bond
				}
				neighbor_type = final_configuration[i];
				//printf("final neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++nsf[neighbor_type - 1]; //number of type A/B/C neighboring atoms in final configuration
					//energy -= nnE[nneA_index[neighbor_type - 1]]*baf[i]; //A-A/B/C bond, told to leave in and comment out
				}
			}
			break;


		case 2:

			for (i=0;i<number_of_neighbors;++i)
			{
				neighbor_type = initial_configuration[i];
				//printf("init neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++nsi[neighbor_type - 1]; //number of type A/B/C neighboring atoms in init configuration
					energy += nnE[nneB_index[neighbor_type - 1]]*baf[i]; //B-A/B/C bond
				}
				neighbor_type = final_configuration[i];
				//printf("final neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++nsf[neighbor_type - 1]; //number of type A/B/C neighboring atoms in final configuration
					//energy -= nnE[nneB_index[neighbor_type - 1]]*baf[i]; //B-A/B/C bond, told to leave in and comment out
				}
			}
			break;

		case 3:
			for (i=0;i<number_of_neighbors;++i)
			{
				neighbor_type = initial_configuration[i];
				//printf("init neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++nsi[neighbor_type - 1]; //number of type A/B/C neighboring atoms in init configuration
					energy += nnE[nneC_index[neighbor_type - 1]]*baf[i]; //C-A/B/C bond
				}
				neighbor_type = final_configuration[i];
				//printf("final neighbor type %d\n", neighbor_type);
				if (neighbor_type > 0) {
					++nsf[neighbor_type - 1]; //number of type A/B/C neighboring atoms in final configuration
					//energy -= nnE[nneC_index[neighbor_type - 1]]*baf[i]; //C-A/B/C bond, told to leave in and comment out
				}
			}
			break;

	}

	//printf("after the sqitch\n");
	for (i=0;i<3;++i)
	{
		ti += nsi[i];
		tf += nsf[i];
	}

	if (ti == 0)
	{
		// this condition corresponds to a diffuser walking through a lattice (a lattice gas)

		energy = -1.0;

	}

	if ((ti > 0)&&(tf <= 1.0))
	{
		energy = 1000.;				// final configuration has no near neighbors,
									// so this effectively corresponds to an evaporation-like event.
									// Don't let it happen!
	}

	// ENHANCE: replace calculating the exp with looking up the value (uhash?) -> speedup?
	//*rate = 1e13*exp(-energy/(kBoltz*temperature)) //+ 1e-4*exp(-(energy-overpotential)/(kBoltz*temperature));
	*rate = 1e13*exp(-energy/(kBoltz*temperature));
	//printf("rate = %le\n", *rate);
	return 0;
}

/******************************************************************************/
/******************************************************************************/

int calculate_evaporation_rate(	int initial_configuration[],			// initial configuration array
									int number_of_neighbors,				// number of 1st near neighbors in current xtal structure
									int atom_type,							// type of atom in consideration for transition
									double nnE[6],							// bond energy (type 2)-(type 2)
									double temperature,						// system temperature
									double overpotential,					// system overpotential
									double *rate)							// return value
{
	int i;
	double energy;

	double EAu = .5;
	double nscale = 0.0;

	static double baf[12] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};			// optional anistropy factor

	energy = 0.0;

	/*if ((ncsk == 1)||(evaporation_flag == false)) //I don't think this happens
	{
		energy = 1000.;
		*rate = 1e4*exp(-(energy)/(kBoltz*temperature));
		return 1;
	}*/

	int nneA_index[3] = {0, 1, 2}; //indices of A-A, A-B, A-C bonds
	int nneB_index[3] = {1, 3, 4}; //indices of B-A, B-B, B-C bonds
	int nneC_index[3] = {2, 4, 5}; //indices of C-A, C-B, C-C bonds
	int neighbor_type; //type of atom for nearest neighbor


	switch(atom_type)
	{
		case 1:
			if (solubility[0]) {
				//A can evaporate
				for (i=0;i<number_of_neighbors;++i)
				{
					neighbor_type = initial_configuration[i];
					if (neighbor_type > 0)
						energy += nnE[nneA_index[neighbor_type - 1]]*baf[i]; //bonding with A
				}
			}
			else
				energy = 1000.; //A cannot evaporate
			*rate = 1e4*exp(-(energy-overpotential)/(kBoltz*temperature));
			break;

		case 2:
			if (solubility[1]) {
				//B can evaporate
				for (i=0;i<number_of_neighbors;++i)
				{
					neighbor_type = initial_configuration[i];
					if (neighbor_type > 0)
						energy += nnE[nneB_index[neighbor_type - 1]]*baf[i]; //bonding with B
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
						energy += nnE[nneC_index[neighbor_type - 1]]*baf[i]; //bonding with C
				}
			}
			else
				energy = 1000.; //C cannot evaporate
			*rate = 1e4*exp(-(energy-overpotential)/(kBoltz*temperature));
			break;
	}
		
	return 0;
}