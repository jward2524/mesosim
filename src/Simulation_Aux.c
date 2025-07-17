#include "Simulation_Aux.h"
#include "Vector.h"
#include "Random.h"
#include "Atoms.h"
#include "Simulation.h"
#include "FileIO.h"
#include "Mesosim.h"
#include <math.h>
#include <errno.h>
#include <stdlib.h>

double ssr;

double normal_x, normal_y, normal_z;

double lhs[6];
double normal_lat[6][3];
int translation_vector[6][3];

/******************************************************************************/
/******************************************************************************/
// XXX: only used for re-deposition
void get_system_rw_radius(struct SimulationEnv* se)
{
	int ss;

	// find minimum axial distance to system edge

	ss = se->system_size_x;
	if (se->system_size_y < ss) ss = se->system_size_y;
	if (se->system_size_z < ss) ss = se->system_size_z;

	ssr = (double)ss/2.;
	ssr = ssr - 5.;
	
	return;
}

/******************************************************************************/
/******************************************************************************/
// updates normal_x, normal_y, normal_z
void get_system_normal(double basis[3][3]) // XXX: supposedly only for vizualization
{
	double nmag;
	double a[3], b[3];
	double n[3];

	// first (x) coordinates of lat vectors
	a[0] = basis[0][0];
	a[1] = basis[1][0];
	a[2] = basis[2][0];

	// second (y) coordinates of lat vectors
	b[0] = basis[0][1];
	b[1] = basis[1][1];
	b[2] = basis[2][1];

	cross(a,b,n);

	normal_x = n[0];
	normal_y = n[1];
	normal_z = n[2];

	nmag = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);

	normal_x = normal_x/nmag;
	normal_y = normal_y/nmag;
	normal_z = normal_z/nmag;

	if (normal_z < 0)
	{
		normal_x *= -1.0;
		normal_y *= -1.0;
		normal_z *= -1.0;
	}

	return;
}

/******************************************************************************/
/******************************************************************************/

//only need to call this when deposition is enabled!
/*int finish_preprocessing(void)
{
	//get_system_rw_radius();
	//get_system_normal();
	return 0;
}*/


/******************************************************************************/
/******************************************************************************/
// updates [iv, iy; primitive_basis, ucell_params, Atoms' cart_coords?; zone_arr, rmat; normal_x, normal_y, normal_z; max_neighbors, se->jump_offset, se->opposite_offset; zi*, zi*shift, *sh], rate_cnt, transition_cnt, atom_cnt, frequency_sum, elapsed_stime, overpotential, next_log_checkpoint
void general_simulation_initialization(struct SimulationState* ss, struct SimulationEnv* se)
{
	// finish initializing structs

	ss->transition_probability.rate_arr_index = (int *)malloc(se->max_rates * sizeof(int));
    ss->transition_probability.lbound = (double *)malloc(se->max_rates * sizeof(double));
    ss->transition_probability.ubound = (double *)malloc(se->max_rates * sizeof(double));

    se->max_transitions = ((MAXIMUM_NUMBER_OF_NEIGHBORS + DISSOLUTION) * se->max_atoms) + 10;
    
    ss->atom_arr = (Atom**) malloc(se->max_atoms * sizeof(Atom*));
    ss->rate_arr = (Rate*) malloc(se->max_rates * sizeof(Rate));
    ss->transition_arr = (Transition**) malloc(se->max_transitions * sizeof(Transition*));

    // null pointer checks
    if (!ss->atom_arr)
    {
        perror("Couldn't allocate memory for atom array");
        clean_and_exit(errno);
    }
    if (!ss->rate_arr)
    {
        perror("Couldn't allocate memory for rate array");
        clean_and_exit(errno);
    }
    if (!ss->transition_arr)
    {
        perror("Couldn't allocate memory for transition array");
        clean_and_exit(errno);
    }

	ss->rate_cnt = 0;	// initialize global transition variables
	ss->transition_cnt = 0;
	ss->atom_cnt = 0;	// initialize global atom variables
	//current_iteration = 0; //not needed if only running 1 simulation at a time // XXX: commented code, never used
	ss->frequency_sum = 0.0;

	ss->elapsed_stime = 0.0;

	ss->overpotential = se->initial_overpotential;

	se->centroid[0] = 0.;
	se->centroid[1] = 0.;
	se->centroid[2] = 0.;

	se->jump_offset = malloc(se->max_neighbors * sizeof(CrystalOffset));
	se->opposite_offset = malloc(se->max_neighbors * sizeof(*se->jump_offset));
	se->ucell_params[0] = 1;
	se->ucell_params[1] = 1;
	se->ucell_params[2] = 1;
	se->ucell_params[3] = 90;
	se->ucell_params[4] = 90;
	se->ucell_params[5] = 90;
	
	if (!se->jump_offset)
    {
        perror("Couldn't allocate memory for jump offset array");
        clean_and_exit(errno);
    }
	if (!se->opposite_offset)
    {
        perror("Couldn't allocate memory for jump offset array");
        clean_and_exit(errno);
    }

	// // XXX: redundant
	// // next_log_checkpoint is initialized to one log_interval_step step
	// if (ls->analysis_type == REGULAR_TIME_INTERVALS) // TODO: reconsider what is happening here
	// 	ls->next_log_checkpoint = ls->next_log_checkpoint;	//do we want this to be true? - overwrites what was in the input file
	// else if (ls->analysis_type == LN_TIME_INTERVALS)
	// ls->next_log_checkpoint = ls->next_log_checkpoint;


	// first, remove any atoms that may exist
	// [ ]: why would atom_cnt not be zero????
	while (ss->atom_cnt != 0) // TODO: start here
		kill_atom(ss->atom_cnt-1, ss, se);

	if (rand_seed > 0) rand_seed = -rand_seed;
	srandj(&rand_seed);
	// atom_cnt=0 for the initialization functions, so some of them end up doing nothing
	get_shifts(se);	// bit shifts for periodic boundary conditions

	// system geometry initialization

	set_primitive_basis(se);
	set_default_orientation(ss->atom_arr, ss->atom_cnt, se->lattice_type, se->rmat, se->primitive_basis); // supposedly was only for visualization
	get_system_normal(se->primitive_basis);	// maybe only for visualization

	// initialize data structures that help figure out which atoms are next to which other atoms

	initialize_neighbor_offsets(se->lattice_type, &se->max_neighbors, se->jump_offset, se->opposite_offset);
	initialize_zones(ss->zone_arr, se);							// initialize zone offsets

	//set_atom_colors(atom_color); // not needed anymore

	return;
}


void initialize_initial_structure(struct SimulationState* ss, struct SimulationEnv* se) // index represents simulation_type, from macros
{
	//printf("I'm in here, simulation index is %d\n", simulation_index);
	switch(se->simulation_type) // TODO: just use the damn macros instead
	{
		case 1:										// flat plane
			initialize_flat_sheet_1(ss, se);
			break;

		case 2:
			initialize_spherical_cluster(ss, se);
			break;
		case 3:
			initialize_from_file(se->atoms_filename); //TODO! THIS IS BIG!
			break;
	}

	// optimizes the atoms added in the initialization routines??
	check_system(ss, se); 
	// organize(ss->atom_arr, ss->atom_cnt, se->primitive_basis);

	return;
}

/********************************************************************************/
/********************************************************************************/
// zi* are the number of zones in that dimension, zi*shift is for bit shifting to find which zone a lattice coordinate corresponds to?
void get_shifts(struct SimulationEnv* se)
{ // updates zi*, zi*shift, *sh
	int temp1;

	temp1 = se->zone_count_u;
	se->zixshift = 0;
	while (temp1 > 1)
	{
		++se->zixshift;
		temp1 = temp1/2;
	}

	temp1 = se->zone_count_v;
	se->ziyshift = 0;
	while (temp1 > 1)
	{
		++se->ziyshift;
		temp1 = temp1/2;
	}

	temp1 = se->zone_count_w;
	se->zizshift = 0;
	while (temp1 > 1)
	{
		++se->zizshift;
		temp1 = temp1/2;
	}

	temp1 = se->system_size_x;
	se->ssxshift = 0;
	while (temp1 > 1)
	{
		++se->ssxshift;
		temp1 = temp1/2;
	}

	temp1 = se->system_size_y;
	se->ssyshift = 0;
	while (temp1 > 1)
	{
		++se->ssyshift;
		temp1 = temp1/2;
	}

	temp1 = se->system_size_z;
	se->sszshift = 0;
	while (temp1 > 1)
	{
		++se->sszshift;
		temp1 = temp1/2;
	}
	// never used, just left and right shift with se->zixshift and se->ssxshift in findzone()
	se->xsh = se->zixshift - se->ssxshift;
	se->ysh = se->ziyshift - se->ssyshift;
	se->zsh = se->zizshift - se->sszshift;
	// TODO: more shifts in zones than in system? se->zone_count_u > se->system_size_x
	return;
}

/******************************************************************************/
/******************************************************************************/

void adjust_pbc(int* u, int* v, int* w, struct SimulationEnv* se) // should be lattice coordinates
{
	// TODO: allow to turn off pbc in a direction
	// if exceeds boundary, either:
	// treat as occupied site - of what composition? atoms will scale box walls like adatoms
	// *treat as not a site - doesn't contribute to energy, can't be transitioned to
	// 		if not a site, return -2 for adjust_pbc, skip atom_at+findzone, don't collect energy or create transition

	// x y z in lattice coordinates
	int u_range = se->simbox_limits_lat[0][1] - se->simbox_limits_lat[0][0];
	int v_range = se->simbox_limits_lat[1][1] - se->simbox_limits_lat[1][0];
	int w_range = se->simbox_limits_lat[2][1] - se->simbox_limits_lat[2][0];

	if (*u < se->simbox_limits_lat[0][0])
		*u += u_range;
	if (*u >= se->simbox_limits_lat[0][1])
		*u -= u_range;

	if (*v < se->simbox_limits_lat[1][0])
		*v += v_range;
	if (*v >= se->simbox_limits_lat[1][1])
		*v -= v_range;

	if (*w < se->simbox_limits_lat[2][0])
		*w += w_range;
	if (*w >= se->simbox_limits_lat[2][1])
		*w -= w_range;

	return;

	// check_pbc(se->primitive_basis)
}

/********************************************************************************/
/********************************************************************************/
// finds the zone indices xy yz zz that correspond to the lattice coordinates xxx yyy zzz
void findzone(int *zone_u, int *zone_v, int *zone_w, int u, int v, int w, struct SimulationEnv* se)
{ 
	// *z are pointers to return indices of the zone, *** are lattice coordinates
	// normalize coordinates to the sblimits, then find which zone
	// (zones / extent) * adjusted_coordinate
	*zone_u = (int) (((double) se->zone_count_u / (se->simbox_limits_lat[0][1] - se->simbox_limits_lat[0][0])) * (u - se->simbox_limits_lat[0][0]));
	*zone_v = (int) (((double) se->zone_count_v / (se->simbox_limits_lat[1][1] - se->simbox_limits_lat[1][0])) * (v - se->simbox_limits_lat[1][0]));
	*zone_w = (int) (((double) se->zone_count_w / (se->simbox_limits_lat[2][1] - se->simbox_limits_lat[2][0])) * (w - se->simbox_limits_lat[2][0]));

	return;
}

/********************************************************************************/
/********************************************************************************/
// updates rmat
void set_default_orientation(Atom** atom_arr, int atom_cnt, int lattice_type, double rmat[3][3], double basis[3][3]) // supposedly for viewing
{
	double axis[3], pnormal[3], axis_mag;
	double zero_point[3] = {0.,0.,0.}, a_a[3];
	double zaxis[3]={0.,0.,1.};
	double zeroa[3]={0.,0.,0.};
	double spin_ax[3];
	double vec_angle;

	// organize(atom_arr, atom_cnt, basis); // atom_cnt ='d 0

	switch(lattice_type)
	{
	    case FCC:
			axis[0] = -1.0;										// {1,1,1}
			axis[1] = 1.0;
			axis[2] = 1.0;
		 	break;

	    case BCC:
			axis[0] = 0.0;										// {1,1,1}
			axis[1] = 1.0;
			axis[2] = -1.0;
	 		break;

	    case SC:
		default:
			axis[0] = 0.0;										// {1,1,1}
			axis[1] = 0.0;
			axis[2] = 1.0;
		 	break;
	}

	axis_mag = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

	// direction cosines of axis normal to plane

	pnormal[0] = axis[0]/axis_mag;
	pnormal[1] = axis[1]/axis_mag;
	pnormal[2] = axis[2]/axis_mag;

	if ((pnormal[0]==zaxis[0])&&(pnormal[1]==zaxis[1])&&(pnormal[1]==zaxis[1]))
	{
		identity2(rmat);
		//removed rotation_notify_flag
		return;
	}

	// orient
	vecdif(pnormal, zero_point, a_a);
	if (((a_a[0]==0.)&&(a_a[1]==0.))||(magnitude(a_a)==0.)) return;

	normto(a_a, zaxis, spin_ax);
	vec_angle = -0.0174533*vangle(zaxis,zeroa,a_a);

	rotmata(spin_ax, vec_angle, rmat);
	transpose(rmat);

	//rotation_notify_flag removed

	return;
}

/********************************************************************************/
/********************************************************************************/
// updates zone (array) based on zi* (zone sizes?), initializes offset to -1 
void initialize_zones(Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], struct SimulationEnv* se)
{
	int i, j, k;
	
	for (i=0;i<se->zone_count_u;++i)
		for (j=0;j<se->zone_count_v;++j)
			for (k=0;k<se->zone_count_w;++k)
				zone_arr[i][j][k].offset = -1;
	return;
}

/********************************************************************************/
/********************************************************************************/
// updates max_neighbors, [se->jump_offset, se->opposite_offset]
void initialize_neighbor_offsets(int lattice_type, int* max_neighbors, CrystalOffset *jump_offset, int *opposite_offset)
{	
	switch(lattice_type)
		{
			case FCC:
				*max_neighbors = 12;
				initialize_jump_offsets(FCC, jump_offset, opposite_offset);
				break;

			case SC:
				*max_neighbors = 6;
				initialize_jump_offsets(SC, jump_offset, opposite_offset);
				break;

			case BCC:
				*max_neighbors = 8;
				initialize_jump_offsets(BCC, jump_offset, opposite_offset);
				break;
		}

	return;
}
// ENHANCE: this is redundant with initialize_neighbor_offsets
// initializes se->jump_offset, se->opposite_offset
void initialize_jump_offsets(int lattice_type, CrystalOffset *jump_offset, int *opposite_offset)	// lattice_type = crystal lattice type
{
	int i;
	int fcc_offs[12] = {11, 10, 7, 4, 3, 6, 5, 2, 9, 8, 1, 0};
	int sc_offs[6] = {1, 0, 3, 2, 5, 4};
	int bcc_offs[8] = {7, 6, 3, 2, 5, 4, 1, 0};

	switch (lattice_type)
	{
		case FCC:
			for (i=0;i<12;++i)
			{
				jump_offset[i].dx = FCC_OFFSET[i].dx;
				jump_offset[i].dy = FCC_OFFSET[i].dy;
				jump_offset[i].dz = FCC_OFFSET[i].dz;
				opposite_offset[i] = fcc_offs[i];
			}

			break;

		case SC:
			for (i=0;i<6;++i)
			{
				jump_offset[i].dx = SC_OFFSET[i].dx;
				jump_offset[i].dy = SC_OFFSET[i].dy;
				jump_offset[i].dz = SC_OFFSET[i].dz;
				opposite_offset[i] = sc_offs[i];
			}

			break;

		case BCC:
			for (i=0;i<8;++i)
			{
				jump_offset[i].dx = BCC_OFFSET[i].dx;
				jump_offset[i].dy = BCC_OFFSET[i].dy;
				jump_offset[i].dz = BCC_OFFSET[i].dz;
				opposite_offset[i] = bcc_offs[i];
			}
			break;
	}

	return;
}

/********************************************************************************/
/********************************************************************************/
// sum of all bond energies in system
void calculate_internal_energy(Atom** atom_arr, int atom_cnt, double* total_internal_energy, struct SimulationEnv* se)
{
	int neighbor, a_type, b_type, bond_idx, env_idx;
	
	*total_internal_energy = 0.;

	// for every atom, for every neighbor, calculate energy
	for (int i = 0; i < atom_cnt; ++i) {
		a_type = atom_arr[i]->type;
		for (int j = 0; j < se->max_neighbors; ++j)
		{
			neighbor = atom_arr[i]->neighbor_atom_idxs[j];

			if (neighbor != -1) //site is not empty
			{
				b_type = atom_arr[neighbor]->type;
				
				//bonds are assumed to be isotropic
				bond_idx = get_bond_index(a_type, b_type, se);
				env_idx = get_env_index(1, bond_idx, se);

				*total_internal_energy += se->nn_energy[env_idx];
			}

		}
	}
	// double counting all bonds, so divide by 2
	*total_internal_energy /= 2.;
	return;
}
/********************************************************************************/
/********************************************************************************/
// updates primitive_basis, ucell_params, Atoms' cart_coords? to match lattice type
void set_primitive_basis(struct SimulationEnv *se) // lattice_type = crystal structure type
{
	switch(se->lattice_type)
	{
		case FCC:
			se->primitive_basis[0][0] = (double)FCCXV1;
			se->primitive_basis[0][1] = (double)FCCXV2;
			se->primitive_basis[0][2] = (double)FCCXV3;

			se->primitive_basis[1][0] = (double)FCCYV1;
			se->primitive_basis[1][1] = (double)FCCYV2;
			se->primitive_basis[1][2] = (double)FCCYV3;

			se->primitive_basis[2][0] = (double)FCCZV1;
			se->primitive_basis[2][1] = (double)FCCZV2;
			se->primitive_basis[2][2] = (double)FCCZV3;
			break;
	
		case BCC:
			se->primitive_basis[0][0] = (double)BCCXV1;
			se->primitive_basis[0][1] = (double)BCCXV2;
			se->primitive_basis[0][2] = (double)BCCXV3;

			se->primitive_basis[1][0] = (double)BCCYV1;
			se->primitive_basis[1][1] = (double)BCCYV2;
			se->primitive_basis[1][2] = (double)BCCYV3;

			se->primitive_basis[2][0] = (double)BCCZV1;
			se->primitive_basis[2][1] = (double)BCCZV2;
			se->primitive_basis[2][2] = (double)BCCZV3;
			break;

		case SC:
			se->primitive_basis[0][0] = (double)SCXV1;
			se->primitive_basis[0][1] = (double)SCXV2;
			se->primitive_basis[0][2] = (double)SCXV3;

			se->primitive_basis[1][0] = (double)SCYV1;
			se->primitive_basis[1][1] = (double)SCYV2;
			se->primitive_basis[1][2] = (double)SCYV3;

			se->primitive_basis[2][0] = (double)SCZV1;
			se->primitive_basis[2][1] = (double)SCZV2;
			se->primitive_basis[2][2] = (double)SCZV3;
			break;
	}

	inver(se->primitive_basis, se->invert_primitive_basis);
	primitive_basis2ucell_params(se->primitive_basis, se->ucell_params);

	// organize(atom_arr, atom_cnt); // ENHANCE: likely unnecessary bc at this point atom_cnt=0
	return;
}

/********************************************************************************/
/********************************************************************************/
// fills initial_config with type of neighbors to atom[at], before jump offset_idx
int get_initial_configuration(int atom_idx, int max_neighbors, Atom** atom_arr, int initial_config[]) // atom_idx is position in atom list, offset_idx is index in se->jump_offset
{	// TODO: rename to remove the 2
   	int offset_idx, neighbor_idx;
	int nn_count = 0; // nearest-neighbors

	for (offset_idx=0; offset_idx<max_neighbors; ++offset_idx)
    {
		neighbor_idx = atom_arr[atom_idx]->neighbor_atom_idxs[offset_idx];

		if (neighbor_idx == -1)
			initial_config[offset_idx] = -1;	// site is empty
	    else
		{
			++nn_count;	// increment number of near neighbors
			initial_config[offset_idx] = atom_arr[neighbor_idx]->type;	// site is occupied by some atom
		}
	}

	return nn_count;
}

/********************************************************************************/
/********************************************************************************/
// fills initial_config with type of neighbors to atom[at], after jump in direction se->jump_offset[offset_idx]
int get_final_configuration(int at, int offset_idx, struct SimulationState *ss, struct SimulationEnv *se, int final_config[]) // offset_idx is position in offset list
{
	int atom_idx;
	int new_x, new_y, new_z;
	int neighbor_x, neighbor_y, neighbor_z;
	
	int nn_cnt = 0; // nearest-neighbors
	// atom position after jump offset_idx
	new_x = ss->atom_arr[at]->lattice[0] + se->jump_offset[offset_idx].dx;
	new_y = ss->atom_arr[at]->lattice[1] + se->jump_offset[offset_idx].dy;
	new_z = ss->atom_arr[at]->lattice[2] + se->jump_offset[offset_idx].dz;

	adjust_pbc(&new_x, &new_y, &new_z, se);

	for (int i=0; i<se->max_neighbors; ++i)
	{
		if (i == se->opposite_offset[offset_idx]) {
			// if direction is where the jump came from, set as empty
			final_config[i] = -1; //hardcode this? // [ ]: is there a case where it won't be empty?
			continue;
		}
		// location of neighbor
		neighbor_x = new_x + se->jump_offset[i].dx;
		neighbor_y = new_y + se->jump_offset[i].dy;
		neighbor_z = new_z + se->jump_offset[i].dz;

		adjust_pbc(&neighbor_x, &neighbor_y, &neighbor_z, se);

		atom_idx = atom_at(neighbor_x, neighbor_y, neighbor_z, ss->atom_arr, ss->zone_arr, se);

		if (atom_idx != -1)
		{
			// if there is an atom present, 'return' its type
			final_config[i] = ss->atom_arr[at]->type;
			++nn_cnt;
		}
		else 
			final_config[i] = -1;
	}

	return nn_cnt;
}

/********************************************************************************/
/********************************************************************************/

void initialize_flat_sheet_1(struct SimulationState *ss, struct SimulationEnv *se)
{
	int i,j,k;
	double nz;
	//printf("Hi there\n");
	for (k = 0; k < se->sheet_thickness; ++k) //new here! loop through z because nothing is buried
	{
		//printf("layer k = %d\n", k);
		for (i=0;i < se->system_size_x; ++i)						// loop through x and y
		{
			for (j=0;j < se->system_size_y; ++j)
			{
				nz = drandj(&rand_seed);
				
				// TODO: make into fxn
				double bar = 0;
				int type = -1;
				do {
					type++;
					bar += se->substrate_composition[type];
				}
				while (bar < nz);
				add_atom(i, j, k, type, NORMAL, ss, se);
			}
		}
	}
	return;
}

/********************************************************************************/
/********************************************************************************/
// ENHANCE: currently adds one extra atom to radius - remove it
void initialize_spherical_cluster(struct SimulationState *ss, struct SimulationEnv *se) // radius of cluster in number of atoms (nearest-neighbor distances)
{
	double center_cart[3]; // cartesian/orthogonal coordinates of center point
	int center_lattice[3]; // lattice coordinates of center point
	int atom_pos_lattice[3]; // atom position in lattice coords
	double atom_pos_cart[3];
	
	double random_num; // random number
	
	double radius_cart; // radius of cluster, in cartesian units

	double dist;

	// for when ss* represents a prism cell
	//center of the cluster is the halfway point
	center_cart[0] = se->system_size_x / 2;
	center_cart[1] = se->system_size_y / 2;
	center_cart[2] = se->system_size_z / 2;
	// TODO: check that the center is at a lattice site?
	cartesian2lattice_site(center_cart, se->invert_primitive_basis, center_lattice);

	// center_lattice[0] = se->system_size_x / 2;
	// center_lattice[1] = se->system_size_y / 2;
	// center_lattice[2] = se->system_size_z / 2;

	// convert lattice distance to cartesian distance using largest (smallest?) lattice vector
	double max_mag = -1; 
	double mag;
	for (int dim = 0; dim < 3; dim++)
	{
		mag = magnitude(se->primitive_basis[dim]);
		if (mag > max_mag)
			max_mag = mag;
	}
	radius_cart = se->cluster_radius * max_mag;

	// algorithm: lattice sphere from cartesian sphere
	// equation: x^2 + y^2 + z^2 <= radius_cart^2
	// convert the 8 corners of the bounding cube into lattice coordinates
	// pick the min and max lattice coordintes from the 6 for each lattice direction
	// loop over lattice coordinates from min to max
	// check if they are in sphere

	// bounding box (bb) limits in cartesian coordinates
	double bblimits_cart[3][2] = {
		{center_cart[0] - radius_cart, center_cart[0] + radius_cart}, // x limits
		{center_cart[1] - radius_cart, center_cart[1] + radius_cart}, // y limits
		{center_cart[2] - radius_cart, center_cart[2] + radius_cart}  // z limits
	};

	// bounding box corners in cartesian coordinates
	double bbcorners_cart[8][3] = {
		{bblimits_cart[0][0], bblimits_cart[1][0], bblimits_cart[2][0]},
		{bblimits_cart[0][0], bblimits_cart[1][0], bblimits_cart[2][1]},
		{bblimits_cart[0][0], bblimits_cart[1][1], bblimits_cart[2][0]},
		{bblimits_cart[0][0], bblimits_cart[1][1], bblimits_cart[2][1]},
		{bblimits_cart[0][1], bblimits_cart[1][0], bblimits_cart[2][0]},
		{bblimits_cart[0][1], bblimits_cart[1][0], bblimits_cart[2][1]},
		{bblimits_cart[0][1], bblimits_cart[1][1], bblimits_cart[2][0]},
		{bblimits_cart[0][1], bblimits_cart[1][1], bblimits_cart[2][1]}
	};

	// convert corners from cartesian coords into atom/lattice coords
	// and find limits in each dimension
	int bblimits_lattice[3][2] = {
		{center_lattice[0], center_lattice[0]}, // u min and max
		{center_lattice[1], center_lattice[1]}, // v min and max
		{center_lattice[2], center_lattice[2]}, // w min and max
	};
	corners2limits(bbcorners_cart, bblimits_lattice, se->invert_primitive_basis);
	
	// check if value exceeds limits for every dimension
	for (int dim_idx = 0; dim_idx < 3; dim_idx++) {
		if (bblimits_cart[dim_idx][0] < 0) {
			printf("ERROR! Spherical cluster passes through periodic boundary conditions\n");
			clean_and_exit(1);
		}
		if (bblimits_cart[dim_idx][1] > (int)(2*center_cart[dim_idx])) {
			printf("ERROR! Spherical cluster passes through periodic boundary conditions\n");
			clean_and_exit(1);
		}
	}

	// iterates through the bounding box of the sphere to identify positions in cluster // ENHANCE: only 52% of loops will be successful - make it more efficient
	int type;
	double bar;
	for (int u = bblimits_lattice[0][0]; u <= bblimits_lattice[0][1]; u++) {
		for (int v = bblimits_lattice[1][0]; v <= bblimits_lattice[1][1]; v++) {
			for (int w = bblimits_lattice[2][0]; w <= bblimits_lattice[2][1]; w++) {
				//convert i, j, k to cartesian coordinates
				atom_pos_lattice[0] = u;
				atom_pos_lattice[1] = v;
				atom_pos_lattice[2] = w;
				
				lattice2cartesian(atom_pos_lattice, se->primitive_basis, atom_pos_cart); // to cartesian coordinates
				dist = 
					(atom_pos_cart[0] - center_cart[0]) * (atom_pos_cart[0] - center_cart[0])
				 	+ (atom_pos_cart[1] - center_cart[1]) * (atom_pos_cart[1] - center_cart[1]) 
					+ (atom_pos_cart[2] - center_cart[2]) * (atom_pos_cart[2] - center_cart[2]); // distance to center

				if (dist <= (radius_cart*radius_cart)) {
					//particle is in bounds
 					random_num = drandj(&rand_seed);
					// determining composition of atom to be placed
					bar = 0;
					type = -1;
					do {
						type++;
						bar += se->substrate_composition[type];
					}
					while (bar < random_num);
					add_atom(u, v, w, type, NORMAL, ss, se);
				}	
			}
		}
	}
	return;
}

/********************************************************************************/
/********************************************************************************/
// get atom positions from a file
void initialize_from_file(char* filename) {
	//does this need more to it?
	// TODO: atom positions from file
	// simulation_parameters_from_file(filename);
	return;
}

// six sides of box
double normal_cart[6][3] = 
{
	{1, 0, 0},
	{0, 1, 0},
	{0, 0, 1},
	{-1, 0, 0},
	{0, -1, 0},
	{0, 0, -1},
};

// int simbox_limits_lat[3][2]; // lattice limits of simulation box in each dimension - for zones
void initialize_simulation_box(struct SimulationEnv* se) //double system_size_x, double system_size_y, double system_size_z)
{
	// assuming simulation box/prism
	// system size in cartesian units [nearest-neighbor (or some other lattice-based) units in cartesian grid]
	// system size of 128 -> 0 to 127, 128th->0
	// need to define box in terms of lattice vectors
	// 6 planes, of form dot(normal, point on plane) = dot(normal, [x,y,z of point to test])
	// normal=(1,0,0); point on plane=(128,0,0) -> 128 = x
	// normal of family (1,0,0) (-1,0,0)
	// point on plane of family  (se->system_size_x,0,0) (0,0,0)
	// x=0 y=0 z=0 x=se->system_size_x y=se->system_size_y z=se->system_size_z
	// to define a region: x>=0 y>=0 z>=0 x<se->system_size_x y<se->system_size_y z<se->system_size_z
	// normals point towards inside of region (keeps the inequality the same)

	// convert into lattice vector form
	// dot(cart2lattice(normal), cart2lattice(point on plane)) {fixed} = dot(cart2lattice(normal), [point to test])
	// [6 planes x 1] vector = [6 normals x 3 coords] @ [3 x 1] = [6 x 1] vector
	// ([6] vector * [6] sign for comparison + [0 or -1] for counting on the boundary) >= 0
	
	// each plane also has associated translation vector, cart2lattice(system size * normal)
	// if outside region, translate by translation vector

	double point_cart[6][3] = 
	{
		{0, 0, 0}, // x
		{0, 0, 0}, // y
		{0, 0, 0}, // z 
		{se->system_size_x, 0, 0}, // x
		{0, se->system_size_y, 0}, // y
		{0, 0, se->system_size_z}, // z
	};
	
	double center_cart[] = {
		se->system_size_x / 2,
		se->system_size_y / 2,
		se->system_size_z / 2,
	};
	int center_lattice[3];
	// TODO: check that the center is at a lattice site
	cartesian2lattice_site(center_cart, se->invert_primitive_basis, center_lattice);

	double sbcorners_cart[8][3] = {
		{0, 0, 0},
		{se->system_size_x, 0, 0},
		{0, se->system_size_y, 0},
		{0, 0, se->system_size_z},
		{se->system_size_x, se->system_size_y, 0},
		{se->system_size_x, 0, se->system_size_z},
		{0, se->system_size_y, se->system_size_z},
		{se->system_size_x, se->system_size_y, se->system_size_z},
	};

	for (int i = 0; i < 3; i++)
	{
		se->simbox_limits_lat[i][0] = center_lattice[0];
		se->simbox_limits_lat[i][1] = center_lattice[0];
	}
	corners2limits(sbcorners_cart, se->simbox_limits_lat, se->invert_primitive_basis);

	// double point_lat[6][3];
	double translation_dist;
	double scaled_normal_cart[3];
	for (int side = 0; side < 6; side++)
	{
		vecmul(normal_cart[side], se->invert_primitive_basis, normal_lat[side]); // doesn't contribute; only for sake of seeing normal in lat coordinates
		// negative normal gives negative dot product
		lhs[side] = fdot(normal_cart[side], point_cart[side]);

		switch (side % 3) {
			case 0:
				translation_dist = se->system_size_x;
				break;
			case 1:
				translation_dist = se->system_size_y;
				break;
			case 2:
				translation_dist = se->system_size_z;
				break;
		}
		// TODO: handle system sizes that don't result in even translation vectors
		// divide scaled_normal_lat by whole number normal_lat (least-common multiple)
		// smallest_translation()
		fconmul(normal_cart[side], translation_dist, scaled_normal_cart, 3);
		cartesian2lattice_site(scaled_normal_cart, se->invert_primitive_basis, translation_vector[side]);
	}
}

// converts cartesian corners to lattice limits along each dimension
// values are rounded towards initial values of limits_lat
void corners2limits(double corners_cart[8][3], int limits_lat[3][2], double inv_basis[3][3])
{
	// convert corners from cartesian coords into atom/lattice coords
	// and find limits in each dimension
	double bbcorners_lattice[8][3];
	// TODO: make data types of cart and lattice coordinates make more sense, and thus the conversion functions
	// find the loop limits from min and max coordinates of corners
	for (int corner_idx = 0; corner_idx < 8; corner_idx++){
		cartesian2lattice(corners_cart[corner_idx], inv_basis, bbcorners_lattice[corner_idx]);

		for (int dim_idx = 0; dim_idx < 3; dim_idx++){
			// if coordinate is smaller than lower limit, change limit (round towards center)
			if (bbcorners_lattice[corner_idx][dim_idx] < (double) limits_lat[dim_idx][0])
				limits_lat[dim_idx][0] = round_towards(bbcorners_lattice[corner_idx][dim_idx], limits_lat[dim_idx][0]);
			
				// if coordinate is largger than upper limit, change limit (round towards center)
			else if (bbcorners_lattice[corner_idx][dim_idx] > (double) limits_lat[dim_idx][1])
				limits_lat[dim_idx][1] = round_towards(bbcorners_lattice[corner_idx][dim_idx], limits_lat[dim_idx][1]);
		}
	}
}

// check if coordinate passed periodic boundary conditions, and translate if it did
void check_pbc(int* u, int* v, int* w, double basis[3][3])
{
	int coords_lat[3] = {*u, *v, *w};
	// double fcoords_lat[3] = {(double) *x, (double) *y, (double) *z};
	int rhs;
	double coords_cart[3];
	lattice2cartesian(coords_lat, basis, coords_cart);
	for (int side = 0; side < 6; side++)
	{
		// rhs contains the point being tested (xyz)
		// negative normal gives negative rhs
		rhs = fdot(normal_cart[side], coords_cart);

		// example conditions to be in region
		// positive normal: 0 (plane) < 2 (coord)
		// negative normal: -128 (plane) < -127 (coord)
		// non-origin sides (3-5) should reroute to 0
		bool cond; 
		if (side >= 3)
			cond = lhs[side] > rhs;
		else
			cond = lhs[side] >= rhs;
		
		if (cond)
		{
			pbc_translate(coords_lat, translation_vector[side]);
			*u = coords_lat[0];
			*v = coords_lat[1];
			*w = coords_lat[2];
		}
	}
}

// translate a point according to boundary conditions
void pbc_translate(int coords_lat[3], int translation_vector[3])
{
	ivecsum(coords_lat, translation_vector, coords_lat);
}

// TODO:
void smallest_translation(double normal_cart[3], double inv_basis[3][3])
{
	// first, get base conversion
	double normal_lat[3];
	vecmul(normal_cart, inv_basis, normal_lat);

	// if conversion has decimals, invert them and try to convert to integers
	// if conversion fails, remove whole number from inversion and try to invert again
	// repeat, multiplying whole numbers together as you go
	// if all integers are the same, multiply normal_lat by integer
	// if not all the same integers, find Least Common Multiple and multipy normal_lat by LCM
}
