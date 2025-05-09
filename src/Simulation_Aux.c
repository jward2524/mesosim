#include "stdafx.h"
#include "Defs.h"
#include "Geometry.h"
#include "Vector.h"
#include "Random.h"
#include "Simulation_Aux.h"
#include "FileIO.h"
#include "Atoms.h"
#include "Simulation.h"

double total_internal_energy = 0;
// [ ]: what are these?
int zixshift, ziyshift, zizshift; // bit shifts for finding zones from coordinates
int ssxshift, ssyshift, sszshift; // used with zi*shift
int zsh, ysh, xsh;	// total bit shifts, zi*shift - ss*shift

// [ ]: what are the units for this? how does it relate to atomic spacing?
int ssx = DSIMSIZE, ssy = DSIMSIZE, ssz = DSIMSIZE;	// system size x, y, z in lattice coordinates
double ssr;
int zix = TTS, ziy = TTS, ziz = TTS;
// defaults are fcc
int lattice_type = FCC;
int number_of_possible_neighbors = 12; // [ ]: this should be dependent on the crystal structure

int sheet_thickness = -1;
int cluster_radius = -1;
char atoms_filename[256] = "";

// [ ]: what is this for?
Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z];

double initialoverpotential = DEFAULT_OVERPOTENTIAL;
double overpotentialramprate = 0.0;
double maxoverpotential = DEFAULT_OVERPOTENTIAL;

double substrate_percent_a = DEFAULT_COMPOSITION_A;
double substrate_percent_b = DEFAULT_COMPOSITION_B;

double initial_logtime = 1.0e-4;

int analysis_type = REGULAR_TIME_INTERVALS;

double logtime_multiplier;

//double vacancy_density = 0.01; //can always add back in

double overpotential_ramp_rate = 0.0;

//int ncsk = 0;

int total_volume_dissolved;

double normal_x, normal_y, normal_z;

/******************************************************************************/
/******************************************************************************/

void get_system_rw_radius(void)
	{
		int ss;

		// find minimum axial distance to system edge

		ss = ssx;
		if (ssy < ss) ss = ssy;
		if (ssz < ss) ss = ssz;

		ssr = (double)ss/2.;
		ssr = ssr - 5.;
		
		return;
	}

/******************************************************************************/
/******************************************************************************/
// updates normal_x, normal_y, normal_z
void get_system_normal(void) // XXX: supposedly only for vizualization
{
	double nmag;
	double a[3], b[3];
	double n[3];

	// first (x) coordinates of lat vectors
	a[0] = primitive_basis[0][0];
	a[1] = primitive_basis[1][0];
	a[2] = primitive_basis[2][0];

	// second (y) coordinates of lat vectors
	b[0] = primitive_basis[0][1];
	b[1] = primitive_basis[1][1];
	b[2] = primitive_basis[2][1];

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
// updates [iv, iy; primitive_basis, ucell_params, Atoms' cart_coords?; rmat; normal_x, normal_y, normal_z; number_of_possible_neighbors, jump_offset, opposite_offset; zi*, zi*shift, *sh], rate_cnt, transition_cnt, atom_cnt, frequency_sum, elapsed_time, overpotential, time_interval_end
void general_simulation_initialization(void)
{
	// first, remove any atoms that may exist
	// [ ]: why would atom_cnt not be zero????
	while (atom_cnt != 0) // TODO: start here
		kill_atom(atom_cnt-1);

	if (rand_seed > 0) rand_seed = -rand_seed;
	srandj(&rand_seed);
	// atom_cnt=0 for the initialization functions, so some of them end up doing nothing
	getshifts();	// bit shifts for periodic boundary conditions

	// system geometry initialization

	set_primitive_basis(lattice_type);
	set_default_orientation(); // supposedly was only for visualization
	get_system_normal();	// maybe only for visualization

	// initialize data structures that help figure out which atoms are next to which other atoms

	initialize_neighbor_offsets();
	initialize_zones();							// initialize zone offsets

	//set_atom_colors(atom_color); // not needed anymore

	rate_cnt = 0;	// initialize global transition variables
	transition_cnt = 0;

	atom_cnt = 0;	// initialize global atom variables
	//current_iteration = 0; //not needed if only running 1 simulation at a time // XXX: commented code, never used
	frequency_sum = 0.0;

	elapsed_time = 0.0;

	overpotential = initialoverpotential;
	// [ ]: how is this related to logging frequency?
	if (analysis_type == REGULAR_TIME_INTERVALS)
		time_interval_end = data_time_interval;	//do we want this to be true?
	else if (analysis_type == LOG_TIME_INTERVALS)
		time_interval_end = initial_logtime;
		
	return;
}


void do_initialize_simulation(int simulation_index) // index represents simulation_type, from macros
{
	//printf("I'm in here, simulation index is %d\n", simulation_index);
	switch(simulation_index) // TODO: just use the damn macros instead
	{
		case 1:										// flat plane
			initialize_flat_sheet_1(sheet_thickness);
			break;

		case 2:
			initialize_spherical_cluster(cluster_radius);
			break;
		case 3:
			initialize_from_file(atoms_filename); //TODO! THIS IS BIG!
			break;
	}
	//printf("My atoms are added\n");
	check_system(); // optimizes the atoms added in the initialization routines
	//printf("My atoms are checked\n");
	organize(atom_arr, atom_cnt);
	//printf("My atoms are organized\n");
	//simulation_initialized = true; //this never really gets used

	return;
}

/********************************************************************************/
/********************************************************************************/
// zi* are the number of zones in that dimension, zi*shift is for bit shifting to find which zone a lattice coordinate corresponds to?
void getshifts(void)
{ // updates zi*, zi*shift, *sh
	int temp1;

	temp1 = zix;
	zixshift = 0;
	while (temp1 > 1)
	{
		++zixshift;
		temp1 = temp1/2;
	}

	temp1 = ziy;
	ziyshift = 0;
	while (temp1 > 1)
	{
		++ziyshift;
		temp1 = temp1/2;
	}

	temp1 = ziz;
	zizshift = 0;
	while (temp1 > 1)
	{
		++zizshift;
		temp1 = temp1/2;
	}

	temp1 = ssx;
	ssxshift = 0;
	while (temp1 > 1)
	{
		++ssxshift;
		temp1 = temp1/2;
	}

	temp1 = ssy;
	ssyshift = 0;
	while (temp1 > 1)
	{
		++ssyshift;
		temp1 = temp1/2;
	}

	temp1 = ssz;
	sszshift = 0;
	while (temp1 > 1)
	{
		++sszshift;
		temp1 = temp1/2;
	}
	// never used, just left and right shift with zixshift and ssxshift in findzone()
	xsh = zixshift - ssxshift;
	ysh = ziyshift - ssyshift;
	zsh = zizshift - sszshift;
	// TODO: more shifts in zones than in system? zix > ssx
	return;
}

/******************************************************************************/
/******************************************************************************/

void adjust_pbc(double *x, double *y, double *z) // lattice coordinates
{
	if (*x < 0.0) *x += (double)ssx;
	if (*x >= (double)ssx) *x -= (double)ssx;

	if (*y < 0.0) *y += (double)ssy;
	if (*y >= (double)ssy) *y -= (double)ssy;

	if (*z < 0.0) *z += (double)ssz;
	if (*z >= (double)ssz) *z -= (double)ssz;

	return;
}

/********************************************************************************/
/********************************************************************************/

void findzone(int *xz, int *yz, int *zz, double xxx, double yyy, double zzz)
{ // *z are pointers to return indices of the zone, *** are lattice coordinates
	int x, y, z;
		
	x = (int)xxx;
	y = (int)yyy;
	z = (int)zzz;
	// TODO: there *may* be undesirable results when using signed integer type and right shifting
	x = x << zixshift;
	x = x >> ssxshift;

	y = y << ziyshift;
	y = y >> ssyshift;

	z = z << zizshift;
	z = z >> sszshift;

	*xz = x;
	*yz = y;
	*zz = z;

	return;
}

/********************************************************************************/
/********************************************************************************/
// updates rmat
void set_default_orientation(void) // supposedly for viewing
{
	static int index[3] = {1,1,1};
	double axis[3], pnormal[3], axis_mag;
	double zero_point[3] = {0.,0.,0.}, a_a[3];
	double zaxis[3]={0.,0.,1.};
	double zeroa[3]={0.,0.,0.};
	double spin_ax[3];
	double vec_angle;

	organize(atom_arr, atom_cnt); // atom_cnt ='d 0

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
void initialize_zones(void)
{
	int i, j, k;
	
	for (i=0;i<zix;++i)
		for (j=0;j<ziy;++j)
			for (k=0;k<ziz;++k)
				zone_arr[i][j][k].offset = -1;
	return;
}

/********************************************************************************/
/********************************************************************************/
// updates number_of_possible_neighbors, [jump_offset, opposite_offset]
void initialize_neighbor_offsets(void)
{	
	switch(lattice_type)
		{
			case FCC:
				number_of_possible_neighbors = 12;
				initialize_jump_offsets(FCC);
				break;

			case SC:
				number_of_possible_neighbors = 6;
				initialize_jump_offsets(SC);
				break;

			case BCC:
				number_of_possible_neighbors = 8;
				initialize_jump_offsets(BCC);
				break;
		}

	return;
}

// initializes jump_offset, opposite_offset
void initialize_jump_offsets(int lattice_type)	// lattice_type = crystal lattice type
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
				jump_offset[i].dx = fcc_offset[i].dx;
				jump_offset[i].dy = fcc_offset[i].dy;
				jump_offset[i].dz = fcc_offset[i].dz;
				opposite_offset[i] = fcc_offs[i];
			}

			break;

		case SC:
			for (i=0;i<6;++i)
			{
				jump_offset[i].dx = sc_offset[i].dx;
				jump_offset[i].dy = sc_offset[i].dy;
				jump_offset[i].dz = sc_offset[i].dz;
				opposite_offset[i] = sc_offs[i];
			}

			break;

		case BCC:
			for (i=0;i<8;++i)
			{
				jump_offset[i].dx = bcc_offset[i].dx;
				jump_offset[i].dy = bcc_offset[i].dy;
				jump_offset[i].dz = bcc_offset[i].dz;
				opposite_offset[i] = bcc_offs[i];
			}
			break;
	}

	return;
}

/********************************************************************************/
/********************************************************************************/
void calculate_internal_energy(int atom_cnt)
{
	int neighbor, type;
	total_internal_energy = 0.;
	int nneA_index[3] = {0, 1, 2}; //indices of A-A, A-B, A-C bonds
	int nneB_index[3] = {1, 3, 4}; //indices of B-A, B-B, B-C bonds
	int nneC_index[3] = {2, 4, 5}; //indices of C-A, C-B, C-C bonds
	for (int i = 0; i < atom_cnt; ++i) {
		for (int j = 0; j < number_of_possible_neighbors; ++j)
		{
			neighbor = atom_arr[i]->occupied_neighbor_sites[j];

			if (neighbor != -1) //site is not empty
			{
				type = atom_arr[neighbor]->type;
				//bonds are assumed to be isotropic
				switch (atom_arr[i]->type) {
					case 1:
						total_internal_energy += nnE[nneA_index[type - 1]];
						break;
					case 2:
						total_internal_energy += nnE[nneB_index[type - 1]];
						break;
					case 3:
						total_internal_energy += nnE[nneC_index[type - 1]];
						break;
					default:
						break;
				}
			}

		}
	}
	total_internal_energy /= 2.;
	return;
}
/********************************************************************************/
/********************************************************************************/
// updates primitive_basis, ucell_params, Atoms' cart_coords? to match lattice type
void set_primitive_basis(int lattice_type) // lattice_type = crystal structure type
{
	switch(lattice_type)
	{
		case FCC:
			primitive_basis[0][0] = (double)FCCXV1;
			primitive_basis[0][1] = (double)FCCXV2;
			primitive_basis[0][2] = (double)FCCXV3;

			primitive_basis[1][0] = (double)FCCYV1;
			primitive_basis[1][1] = (double)FCCYV2;
			primitive_basis[1][2] = (double)FCCYV3;

			primitive_basis[2][0] = (double)FCCZV1;
			primitive_basis[2][1] = (double)FCCZV2;
			primitive_basis[2][2] = (double)FCCZV3;
			break;
	
		case BCC:
			primitive_basis[0][0] = (double)BCCXV1;
			primitive_basis[0][1] = (double)BCCXV2;
			primitive_basis[0][2] = (double)BCCXV3;

			primitive_basis[1][0] = (double)BCCYV1;
			primitive_basis[1][1] = (double)BCCYV2;
			primitive_basis[1][2] = (double)BCCYV3;

			primitive_basis[2][0] = (double)BCCZV1;
			primitive_basis[2][1] = (double)BCCZV2;
			primitive_basis[2][2] = (double)BCCZV3;
			break;

		case SC:
			primitive_basis[0][0] = (double)SCXV1;
			primitive_basis[0][1] = (double)SCXV2;
			primitive_basis[0][2] = (double)SCXV3;

			primitive_basis[1][0] = (double)SCYV1;
			primitive_basis[1][1] = (double)SCYV2;
			primitive_basis[1][2] = (double)SCYV3;

			primitive_basis[2][0] = (double)SCZV1;
			primitive_basis[2][1] = (double)SCZV2;
			primitive_basis[2][2] = (double)SCZV3;
			break;
	}

	inver(primitive_basis, invert_primitive_basis);
	primitive_basis2ucell_params(primitive_basis, ucell_params);

	organize(atom_arr, atom_cnt); // ENHACNE: likely unnecessary bc at this point atom_cnt=0
	return;
}

/********************************************************************************/
/********************************************************************************/
// fills initial_config with type of neighbors to atom[at], before jump offset_idx
int get_initial_configuration2(int atom_idx, int offset_idx, int initial_config[]) // atom_idx is position in atom list, offset_idx is index in jump_offset
{	// TODO: rename to remove the 2
   	int i, j;
	int nn_count = 0; // nearest-neighbors

	for (i=0; i<number_of_possible_neighbors; ++i)
    {
		j = atom_arr[atom_idx]->occupied_neighbor_sites[i];

		if (j == -1)
			initial_config[i] = -1;	// site is empty
	    else
		{
			++nn_count;	// increment number of near neighbors
			initial_config[i] = atom_arr[atom_idx]->type;	// site is occupied by some atom
		}
	}

	return nn_count;
}

/********************************************************************************/
/********************************************************************************/
// fills initial_config with type of neighbors to atom[at], after jump in direction jump_offset[offset_idx]
int get_final_configuration2(int at, int offset_idx, int final_config[]) // offset_idx is position in offset list
{
	int i, j, k;
	double new_x, new_y, new_z;
	double neighbor_x, neighbor_y, neighbor_z;
	//int n = 0; // XXX:
	int nn_cnt = 0; // nearest-neighbors
	// atom position after jump offset_idx
	new_x = atom_arr[at]->lattice[0] + jump_offset[offset_idx].dx;
	new_y = atom_arr[at]->lattice[1] + jump_offset[offset_idx].dy;
	new_z = atom_arr[at]->lattice[2] + jump_offset[offset_idx].dz;

	//printf("before pbc xyz %lf %lf %lf\n", x, y, z); // XXX: commented print
	adjust_pbc(&new_x, &new_y, &new_z);

	//printf("after pbc xyz %lf %lf %lf\n", x, y, z);
	for (i=0; i<number_of_possible_neighbors; ++i)
      	{
			//printf("offset_idx = %d, i = %d\n", offset_idx, i);
			if (i == opposite_offset[offset_idx]) { // if direction is where the jump came from, set as empty 
				final_config[i] = -1; //hardcode this? // [ ]: is there a case where it won't be empty?
				//printf("opposite offset! final_config[i] = %d\n", final_config[i]);
				continue;
			}
			// location of neighbor
      		neighbor_x = new_x + jump_offset[i].dx;
			neighbor_y = new_y + jump_offset[i].dy;
	        neighbor_z = new_z + jump_offset[i].dz;

			//printf("before pbc nxyz %lf %lf %lf\n", neighbor_x, neighbor_y, neighbor_z);
			adjust_pbc(&neighbor_x, &neighbor_y, &neighbor_z);
			//printf("after pbc nxyz %lf %lf %lf\n", neighbor_x, neighbor_y, neighbor_z);
	        j = atom_at(neighbor_x, neighbor_y, neighbor_z);
			//printf("j = %d\n", j);
	        if (j != -1)
		        { // if there is an atom present, 'return' its type
					final_config[i] = atom_arr[at]->type;
					//printf("at = %d, atom[at]->type = %d, atom[j]->type = %d\n", at, atom[at]->type, atom[j]->type);
					++nn_cnt;
				}
			else final_config[i] = -1;
		}
	/*printf("final config: "); // XXX: commented print
	for (i=0;i<number_of_possible_neighbors;++i)
		printf("%d ", final_config[i]);
	printf("\n");*/

	return nn_cnt;
}

/********************************************************************************/
/********************************************************************************/

void initialize_flat_sheet_1(int z)
{
	int i,j,k;
	double nz;
	//printf("Hi there\n");
	for (k = 0; k < z; ++k) //new here! loop through z because nothing is buried
	{
		//printf("layer k = %d\n", k);
		for (i=0;i<ssx;++i)						// loop through x and y
		{
			for (j=0;j<ssy;++j)
			{
				nz = drandj(&rand_seed);
				//printf("i, j, k = %d, %d, %d\n", i, j, k);
				if (nz <= substrate_percent_a)
					add_atom(i, j, k, 1, NORMAL);
				else if (nz <= substrate_percent_a + substrate_percent_b)
					add_atom(i, j, k, 2, NORMAL);
				else
					add_atom(i, j, k, 3, NORMAL);
			}
		}
	}
	return;
}

/********************************************************************************/
/********************************************************************************/

void initialize_spherical_cluster(int radius_lattice) // radius of cluster in number of atoms
{
	double center_cart[3]; // cartesian/orthogonal coordinates of center point
	double center_lattice[3]; // lattice coordinates of center point
	double atom_pos_lattice[3], atom_pos_cart[3]; // atom position in lattice coords (atom_pos_lattice), cartesian/orthogonal coords (atom_pos_cart)
	
	double random_num; // random number
	
	double radius_cart; // radius of cluster, in cartesian units

	double dist;

	//center of the cluster is the halfway point - center_lattice=lattice/atom coordinate of cluster center
	center_lattice[0] = ssx/2.;
	center_lattice[1] = ssy/2.;
	center_lattice[2] = ssz/2.;
	// # of atoms * translation vector -> cartesian
	vecmul(center_lattice, primitive_basis, center_cart);	// center_cart is the cartesian coordinates of the central point
	// BUG: center_cart is twice what I think it should be - primitive_basis isn't normalized (not unit vectors)
	radius_cart = (double)radius_lattice; // BUG: no shot this is right; max_lattice - center_lattice neq radius_lattice; should be radius_lattice * mag([largest?] lattice vector)

	// lattice sphere from cartesian sphere - algorithm
	// equation: x^2 + y^2 + z^2 <= radius_cart^2
	// convert the 8 corners of the bounding cube into lattice coordinates
	// pick the min and max lattice coordintes from the 6 for each lattice direction
	// loop over lattice coordinates from min to max
	// check if they are in sphere

	// bounding box (bb) limits in cartesian coordinates
	int bblimits_cart[3][2] = {
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
	double bbcorners_lattice[8][3];
	double bblimits_lattice[3][2] = {
		{0.0, 0.0},
		{0.0, 0.0},
		{0.0, 0.0}
	};
	// TODO: make data types of cart and lattice coordinates make more sense, and thus the conversion functions
	for (int corner_idx = 0; corner_idx < 8; corner_idx++){
		vecmul(bbcorners_cart[corner_idx], invert_primitive_basis, bbcorners_lattice[corner_idx]);

		// check if value exceeds limits for every dimension
		for (int dim_idx = 0; dim_idx < 3; dim_idx++){
			if (bbcorners_lattice[corner_idx][dim_idx] < bblimits_lattice[dim_idx][0])
				bblimits_lattice[dim_idx][0] = bbcorners_lattice[corner_idx][dim_idx];
			else if (bbcorners_lattice[corner_idx][dim_idx] > bblimits_lattice[dim_idx][1])
				bblimits_lattice[dim_idx][1] = bbcorners_lattice[corner_idx][dim_idx];
		}
	}

	for (int dim_idx = 0; dim_idx < 3; dim_idx++) {
		if (bblimits_lattice[dim_idx][0] < 0) {
			printf("ERROR! Spherical cluster passes through periodic boundary conditions\n");
			return;
		}
		if ((int)bblimits_lattice[dim_idx][1] > (int)(2*center_lattice[dim_idx])) {
			printf("ERROR! Spherical cluster passes through periodic boundary conditions\n");
			return;
		}
	}

	// iterates through the bounding box of the sphere to identify positions in cluster // ENHANCE: only 52% of loops will be successful - make it more efficient
	for (int u = bblimits_lattice[0][0]; u <= bblimits_lattice[0][1]; u++) {
		for (int v = bblimits_lattice[1][0]; v <= bblimits_lattice[1][1]; v++) {
			for (int w = bblimits_lattice[2][0]; w <= bblimits_lattice[2][1]; w++) {
				//convert i, j, k to cartesian coordinates
				atom_pos_lattice[0] = u;
				atom_pos_lattice[1] = v;
				atom_pos_lattice[2] = w;
				
				vecmul(atom_pos_lattice, primitive_basis, atom_pos_cart); // to cartesian coordinates
				dist = 
					(atom_pos_cart[0] - center_cart[0]) * (atom_pos_cart[0] - center_cart[0])
				 	+ (atom_pos_cart[1] - center_cart[1]) * (atom_pos_cart[1] - center_cart[1]) 
					+ (atom_pos_cart[2] - center_cart[2]) * (atom_pos_cart[2] - center_cart[2]); // distance to center

				if (dist <= (radius_cart*radius_cart)) {
					//particle is in bounds
 					random_num = drandj(&rand_seed);
					// determining composition of atom to be placed
					if (random_num < substrate_percent_a)
						add_atom(u, v, w, 1, NORMAL);
					else if (random_num < substrate_percent_a + substrate_percent_b)
						add_atom(u, v, w, 2, NORMAL);
					else
						add_atom(u, v, w, 3, NORMAL);
				}	
			}
		}
	}
	// conversion from aotm->lattice to cartesian and store in atom->cart_coord
	organize(atom_arr, atom_cnt);

	return;
}

/********************************************************************************/
/********************************************************************************/

void initialize_from_file(char* filename) {
	//does this need more to it?
	get_input_file(filename);
	return;
}