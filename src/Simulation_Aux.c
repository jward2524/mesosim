#include "stdafx.h"
//#include "Mesosim Resources.h"
#include "Defs.h"
#include "Geometry.h"
#include "Prototypes.h"
#include "Vector.h"
#include "Global_Externs.h"
#include "Simulation_Global_Externs.h"

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

void get_system_normal(void)
{
	double nmag;
	double a[3], b[3];
	double n[3];

	a[0] = latmat[0][0];
	a[1] = latmat[1][0];
	a[2] = latmat[2][0];

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

	while (nat != 0)
		kill_atom(nat-1);

	if (seed > 0) seed = -seed;
	srandj(&seed);

	getshifts();								// bit shifts for periodic boundary conditions

	// system geometry initialization

	set_latmat(lattice_type);
	set_default_orientation();
	get_system_normal();

	// initialize data structures that help figure out which atoms are next to which other atoms

	initialize_neighbor_offsets();
	initialize_zones();							// initialize zone offsets

	//set_atom_colors(atom_color); // not needed anymore

	number_rates = 0;							// initialize global transition variables
	total_current_transitions = 0;

	nat = 0;									// initialize global atom variables
	//current_iteration = 0; //not needed if only running 1 simulation at a time
	sum_of_frequencies = 0.0;

	elapsed_time = 0.0;

	overpotential = initialoverpotential;

	if (analysis_type == REGULAR_TIME_INTERVALS)
		time_interval_end = data_time_interval; //do we want this to be true?
	else if (analysis_type == LOG_TIME_INTERVALS)
		time_interval_end = initial_logtime;	
		
	return;
}


void do_initialize_simulation(int simulation_index)
{
	int i, j, n, m;
	//printf("I'm in here, simulation index is %d\n", simulation_index);
	switch(simulation_index)
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
	check_system();					// optimizes the atoms added in the initialization routines
	//printf("My atoms are checked\n");
	organize(atom, nat);
	//printf("My atoms are organized\n");
	//simulation_initialized = true; //this never really gets used

	return;
}

/********************************************************************************/
/********************************************************************************/

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

	xsh = zixshift - ssxshift;
	ysh = ziyshift - ssyshift;
	zsh = zizshift - sszshift;
	
	return;
}

/******************************************************************************/
/******************************************************************************/

void adjust_pbc(double *x, double *y, double *z)
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
{
	int x, y, z;
		
	x = (int)xxx;
	y = (int)yyy;
	z = (int)zzz;

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

void set_default_orientation(void)
{
	static int index[3] = {1,1,1};
	double axis[3], pnormal[3], axismag;
	double zeropoint[3] = {0.,0.,0.}, a_a[3];
	double zaxis[3]={0.,0.,1.};
	double zeroa[3]={0.,0.,0.};
	double spinax[3];
	double vanle;

	organize(atom, nat);

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


void initialize_jump_offsets(int l_t)
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

void set_latmat(int lt)
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

	organize(atom, nat);
	return;
}

/********************************************************************************/
/********************************************************************************/

int get_initial_configuration2(int at, int vc, int initial_config[])
{
   	int i, j;
	int nn = 0;

	for (i=0;i<number_of_possible_neighbors;++i)
    {
		j = atom[at]->occupied_neighbor_sites[i];

		if (j == -1)
			initial_config[i] = -1;							// site is empty
	    else
		{
			++nn;										// increment number of near neighbors
			initial_config[i] = atom[at]->type;	// site is occupied by some atom
		}
	}

	return nn;
}

/********************************************************************************/
/********************************************************************************/

int get_final_configuration2(int at, int vc, int final_config[])
{
	int i, j, k;
	double x, y, z;
	double nx, ny, nz;
	//int n = 0;
	int nnn = 0;

	x = atom[at]->lattice[0] + jump_offset[vc].dx;
	y = atom[at]->lattice[1] + jump_offset[vc].dy;
	z = atom[at]->lattice[2] + jump_offset[vc].dz;

	//printf("before pbc xyz %lf %lf %lf\n", x, y, z);
	adjust_pbc(&x, &y, &z);

	//printf("after pbc xyz %lf %lf %lf\n", x, y, z);
	for (i=0;i<number_of_possible_neighbors;++i)
      	{
			//printf("vc = %d, i = %d\n", vc, i);
			if (i == opposite_offset[vc]) {
				final_config[i] = -1; //hardcode this?
				//printf("opposite offset! final_config[i] = %d\n", final_config[i]);
				continue;
			}

      		nx = x + jump_offset[i].dx;
			ny = y + jump_offset[i].dy;
	        nz = z + jump_offset[i].dz;

			//printf("before pbc nxyz %lf %lf %lf\n", nx, ny, nz);
			adjust_pbc(&nx, &ny, &nz);
			//printf("after pbc nxyz %lf %lf %lf\n", nx, ny, nz);
	        j = atom_at(nx, ny, nz);
			//printf("j = %d\n", j);
	        if (j != -1)
		        {
					final_config[i] = atom[at]->type;
					//printf("at = %d, atom[at]->type = %d, atom[j]->type = %d\n", at, atom[at]->type, atom[j]->type);
					++nnn;
				}
			else final_config[i] = -1;
		}
	/*printf("final config: ");
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

void initialize_spherical_cluster_1(int radius_of_sphere)
{
	int i,j,k;

	double c[3];
	double lc[3];
	double p[3], op[3];
	
	double rn;
	
	double radius;

	int sa;

	double min_xyz[3], max_xyz[3];
	double min_lat[3], max_lat[3];
	double dist;

	//center of the cluster is the halfway point
	lc[0] = ssx/2.;
	lc[1] = ssy/2.;
	lc[2] = ssz/2.;

	vecmul(lc, latmat, c);				// c is the orthogonal coordinates of the central point

	radius = (double)radius_of_sphere;

	//find the min/max x, y, z points
	for (i = 0; i < 3; ++i) {
		min_xyz[i] = c[i] - radius;
		max_xyz[i] = c[i] + radius;
	}

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
				vecmul(p, latmat, op);
				dist = (op[0] - c[0]) * (op[0] - c[0]) + (op[1] - c[1]) * (op[1] - c[1]) + (op[2] - c[2]) * (op[2] - c[2]);
				if (dist <= (radius*radius)) {
					//particle is in bounds
					rn = drandj(&seed);

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