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
Zone zone[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z];

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

void get_system_normal(void) // XXX: supposedly only for vizualization
{
	double nmag;
	double a[3], b[3];
	double n[3];

	// first (x) coordinates of lat vectors
	a[0] = latmat[0][0];
	a[1] = latmat[1][0];
	a[2] = latmat[2][0];

	// second (y) coordinates of lat vectors
	b[0] = latmat[0][1];
	b[1] = latmat[1][1];
	b[2] = latmat[2][1];

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

void general_simulation_initialization(void)
{
	int i;

	// first, remove any atoms that may exist
	// [ ]: why would nat not be zero????
	while (nat != 0)
		kill_atom(nat-1);

	if (seed > 0) seed = -seed;
	srandj(&seed);
	// nat=0 for the initialization functions, so some of them end up doing nothing
	getshifts();	// bit shifts for periodic boundary conditions

	// system geometry initialization

	set_latmat(lattice_type);
	set_default_orientation(); // supposedly was only for visualization
	get_system_normal();	// maybe only for visualization

	// initialize data structures that help figure out which atoms are next to which other atoms

	initialize_neighbor_offsets();
	initialize_zones();							// initialize zone offsets

	//set_atom_colors(atom_color); // not needed anymore

	number_rates = 0;							// initialize global transition variables
	total_current_transitions = 0;

	nat = 0;									// initialize global atom variables
	//current_iteration = 0; //not needed if only running 1 simulation at a time // XXX: commented code, never used
	sum_of_frequencies = 0.0;

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
	int i, j, n, m;
	//printf("I'm in here, simulation index is %d\n", simulation_index);
	switch(simulation_index) // TODO: just use the damn macros instead
	{
		case 1:										// flat plane
			initialize_flat_sheet_1(sheet_thickness);
			break;

		case 2:
			initialize_spherical_cluster_1(cluster_radius);
			break;
		case 3:
			initialize_from_file(atoms_filename); //TODO! THIS IS BIG!
			break;
	}
	//printf("My atoms are added\n");
	check_system(); // optimizes the atoms added in the initialization routines
	//printf("My atoms are checked\n");
	organize(atom, nat);
	//printf("My atoms are organized\n");
	//simulation_initialized = true; //this never really gets used

	return;
}

/********************************************************************************/
/********************************************************************************/
// zi* are the number of zones in that dimension, zi*shift is for bit shifting to find which zone a lattice coordinate corresponds to?
void getshifts(void)
{
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

void set_default_orientation(void) // supposedly for viewing
{
	static int index[3] = {1,1,1};
	double axis[3], pnormal[3], axismag;
	double zeropoint[3] = {0.,0.,0.}, a_a[3];
	double zaxis[3]={0.,0.,1.};
	double zeroa[3]={0.,0.,0.};
	double spinax[3];
	double vanle;

	organize(atom, nat); // nat ='d 0

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

	axismag = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

	// direction cosines of axis normal to plane

	pnormal[0] = axis[0]/axismag;
	pnormal[1] = axis[1]/axismag;
	pnormal[2] = axis[2]/axismag;

	if ((pnormal[0]==zaxis[0])&&(pnormal[1]==zaxis[1])&&(pnormal[1]==zaxis[1]))
	{
		identity2(rmat);
		//removed rotation_notify_flag
		return;
	}

	// orient
	vecdif(pnormal, zeropoint, a_a);
	if (((a_a[0]==0.)&&(a_a[1]==0.))||(magnitude(a_a)==0.)) return;

	normto(a_a, zaxis, spinax);
	vanle = -0.0174533*vangle(zaxis,zeroa,a_a);

	rotmata(spinax,vanle,rmat);
	transpose(rmat);

	//rotation_notify_flag removed

	return;
}

/********************************************************************************/
/********************************************************************************/
// creates zone array based on zi* (zone sizes?), initializes offset to -1; 
void initialize_zones(void)
{
	int i, j, k;
	
	for (i=0;i<zix;++i)
		for (j=0;j<ziy;++j)
			for (k=0;k<ziz;++k)
				zone[i][j][k].offset = -1;
	return;
}

/********************************************************************************/
/********************************************************************************/

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


void initialize_jump_offsets(int l_t)	// l_t = lattice type
	{
		int i;
		int fcc_offs[12] = {11, 10, 7, 4, 3, 6, 5, 2, 9, 8, 1, 0};
		int sc_offs[6] = {1, 0, 3, 2, 5, 4};
		int bcc_offs[8] = {7, 6, 3, 2, 5, 4, 1, 0};

		switch (l_t)
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
void calculate_internal_energy(int nat) {
	int neighbor, type;
	total_internal_energy = 0.;
	int nneA_index[3] = {0, 1, 2}; //indices of A-A, A-B, A-C bonds
	int nneB_index[3] = {1, 3, 4}; //indices of B-A, B-B, B-C bonds
	int nneC_index[3] = {2, 4, 5}; //indices of C-A, C-B, C-C bonds
	for (int i = 0; i < nat; ++i) {
		for (int j = 0; j < number_of_possible_neighbors; ++j)
		{
			neighbor = atom[i]->occupied_neighbor_sites[j];

			if (neighbor != -1) //site is not empty
			{
				type = atom[neighbor]->type;
				//bonds are assumed to be isotropic
				switch (atom[i]->type) {
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

void set_latmat(int lt) // lt = lattice type
{
	switch(lt)
	{
		case FCC:
			latmat[0][0] = (double)FCCXV1;
			latmat[0][1] = (double)FCCXV2;
			latmat[0][2] = (double)FCCXV3;

			latmat[1][0] = (double)FCCYV1;
			latmat[1][1] = (double)FCCYV2;
			latmat[1][2] = (double)FCCYV3;

			latmat[2][0] = (double)FCCZV1;
			latmat[2][1] = (double)FCCZV2;
			latmat[2][2] = (double)FCCZV3;
			break;
	
		case BCC:
			latmat[0][0] = (double)BCCXV1;
			latmat[0][1] = (double)BCCXV2;
			latmat[0][2] = (double)BCCXV3;

			latmat[1][0] = (double)BCCYV1;
			latmat[1][1] = (double)BCCYV2;
			latmat[1][2] = (double)BCCYV3;

			latmat[2][0] = (double)BCCZV1;
			latmat[2][1] = (double)BCCZV2;
			latmat[2][2] = (double)BCCZV3;
			break;

		case SC:
			latmat[0][0] = (double)SCXV1;
			latmat[0][1] = (double)SCXV2;
			latmat[0][2] = (double)SCXV3;

			latmat[1][0] = (double)SCYV1;
			latmat[1][1] = (double)SCYV2;
			latmat[1][2] = (double)SCYV3;

			latmat[2][0] = (double)SCZV1;
			latmat[2][1] = (double)SCZV2;
			latmat[2][2] = (double)SCZV3;
			break;
	}

	inver(latmat, ilatmat);
	latmat_to_cell(latmat, cell);

	organize(atom, nat); // ENHACNE: likely unnecessary bc at this point nat=0
	return;
}

/********************************************************************************/
/********************************************************************************/
// fills initial_config with type of neighbors to atom[at], before jump vc
int get_initial_configuration2(int at, int vc, int initial_config[]) // at is position in atom list, vc is index in offset list?
{	// TODO: rename to remove the 2
   	int i, j;
	int nn = 0;

	for (i=0;i<number_of_possible_neighbors;++i)
    {
		j = atom[at]->occupied_neighbor_sites[i];

		if (j == -1)
			initial_config[i] = -1;	// site is empty
	    else
		{
			++nn;	// increment number of near neighbors
			initial_config[i] = atom[at]->type;	// site is occupied by some atom
		}
	}

	return nn;
}

/********************************************************************************/
/********************************************************************************/
// fills initial_config with type of neighbors to atom[at], after jump in direction jump_offset[vc]
int get_final_configuration2(int at, int vc, int final_config[]) // vc is position in offset list
{
	int i, j, k;
	double x, y, z;
	double nx, ny, nz;
	//int n = 0;
	int nnn = 0;
	// atom position after jump vc
	x = atom[at]->lattice[0] + jump_offset[vc].dx;
	y = atom[at]->lattice[1] + jump_offset[vc].dy;
	z = atom[at]->lattice[2] + jump_offset[vc].dz;

	//printf("before pbc xyz %lf %lf %lf\n", x, y, z); // XXX: commented print
	adjust_pbc(&x, &y, &z);

	//printf("after pbc xyz %lf %lf %lf\n", x, y, z);
	for (i=0;i<number_of_possible_neighbors;++i)
      	{
			//printf("vc = %d, i = %d\n", vc, i);
			if (i == opposite_offset[vc]) { // if direction is where the jump came from, set as empty 
				final_config[i] = -1; //hardcode this? // [ ]: is there a case where it won't be empty?
				//printf("opposite offset! final_config[i] = %d\n", final_config[i]);
				continue;
			}
			// location of neighbor
      		nx = x + jump_offset[i].dx;
			ny = y + jump_offset[i].dy;
	        nz = z + jump_offset[i].dz;

			//printf("before pbc nxyz %lf %lf %lf\n", nx, ny, nz);
			adjust_pbc(&nx, &ny, &nz);
			//printf("after pbc nxyz %lf %lf %lf\n", nx, ny, nz);
	        j = atom_at(nx, ny, nz);
			//printf("j = %d\n", j);
	        if (j != -1)
		        { // if there is an atom present, 'return' its type
					final_config[i] = atom[at]->type;
					//printf("at = %d, atom[at]->type = %d, atom[j]->type = %d\n", at, atom[at]->type, atom[j]->type);
					++nnn;
				}
			else final_config[i] = -1;
		}
	/*printf("final config: "); // XXX: commented print
	for (i=0;i<number_of_possible_neighbors;++i)
		printf("%d ", final_config[i]);
	printf("\n");*/

	return nnn;
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
				nz = drandj(&seed);
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

void initialize_spherical_cluster_1(int radius_of_sphere) // radius of cluster in number of atoms
{
	int i,j,k;

	double c[3]; // cartesian/orthogonal coordinates of center point
	double lc[3]; // lattice coordinates of center point
	double p[3], op[3]; // atom position in lattice coords (p), cartesian/orthogonal coords (op)
	
	double rn; // random number
	
	double radius; // radius of cluster, in cartesian units

	int sa; // substrate atom? the atom created // XXX: unused

	double min_xyz[3], max_xyz[3];
	double min_lat[3], max_lat[3];
	double dist;

	//center of the cluster is the halfway point - lc=lattice/atom coordinate of cluster center
	lc[0] = ssx/2.;
	lc[1] = ssy/2.;
	lc[2] = ssz/2.;
	// # of atoms * translation vector
	vecmul(lc, latmat, c);	// c is the cartesian coordinates of the central point
	// BUG: c is twice what I think it should be - latmat isn't normalized (not unit vectors)
	radius = (double)radius_of_sphere;
	// TODO: change from min/max in xyz to equation of a sphere; radial coordinates?
	// find the min/max x, y, z points (cartesian coords)
	for (i = 0; i < 3; ++i) {
		min_xyz[i] = c[i] - radius;
		max_xyz[i] = c[i] + radius;
	}
	// turn min/max from cartesian coords into atom/lattice coords
	vecmul(min_xyz, ilatmat, min_lat);
	vecmul(max_xyz, ilatmat, max_lat);

	for (i = 0; i < 3; ++i) {
		min_lat[i] = (int)min_lat[i];
		if (min_lat[i] < 0) {
			printf("ERROR! Spherical cluster passes through periodic boundary conditions\n");
			return;
		}
		max_lat[i] = (int)max_lat[i];
		if (max_lat[i] > (int)(2*lc[i])) {
			printf("ERROR! Spherical cluster passes through periodic boundary conditions\n");
			return;
		}
	}

	for (i = min_lat[0]; i <= max_lat[0]; ++i) {
		for (j = min_lat[1]; j <= max_lat[1]; ++j) {
			for (k = min_lat[2]; k <= max_lat[2]; ++k) {
				//convert i, j, k to cartesian coordinates
				p[0] = i;
				p[1] = j;
				p[2] = k;
				vecmul(p, latmat, op); // to cartesian coordinates
				dist = (op[0] - c[0]) * (op[0] - c[0]) + (op[1] - c[1]) * (op[1] - c[1]) + (op[2] - c[2]) * (op[2] - c[2]); // distance to center
				if (dist <= (radius*radius)) {
					//particle is in bounds
					rn = drandj(&seed);
					// determining composition of atom to be placed
					if (rn < substrate_percent_a)
						sa = add_atom(i, j, k, 1, NORMAL);
					else if (rn < substrate_percent_a + substrate_percent_b)
						sa = add_atom(i, j, k, 2, NORMAL);
					else
						sa = add_atom(i, j, k, 3, NORMAL);
				}	
			}
		}
	}

	organize(atom, nat);

	return;
}

/********************************************************************************/
/********************************************************************************/

void initialize_from_file(char* filename) {
	//does this need more to it?
	get_input_file(filename);
	return;
}