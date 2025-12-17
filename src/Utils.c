#include "Utils.h"
#include "Vector.h"
#include "math.h"
#include <assert.h>

double ssr;

// for simulation box
static double lhs[6];
static int translation_vector[6][3];

// six sides of box
double normal_cart[6][3] = 
{
	{1., 0., 0.},
	{0., 1., 0.},
	{0., 0., 1.},
	{-1., 0., 0.},
	{0., -1., 0.},
	{0., 0., -1.},
};

const double DEFAULT_EPSILON = 1e-3;

static int int_check(double fvalue, int ireference, double epsilon);

// checks if double is within range of an integer - returns 1 for true, 0 for false
static int int_check(double fvalue, int ireference, double epsilon){
	return fabs(fvalue - (double) ireference) < epsilon;
}

void lattice2int(double fcoords[3], int coords[3], double epsilon){
	for (int dim_idx = 0; dim_idx < 3; dim_idx++){
		double u = fcoords[dim_idx];
		int comp = round(u);
		int res = int_check(u, comp, epsilon);
		assert(res == 1);
		coords[dim_idx] = comp;
	}
}

// convert a site from cartesian coordinates to lattice coordinates
void cartesian2lattice_site(double ccart[3], double invert_primitive_basis[3][3], int clattice[3]){
	double fclattice[3];
	vecmul(ccart, invert_primitive_basis, fclattice);
	lattice2int(fclattice, clattice, DEFAULT_EPSILON);
}

void cartesian2lattice(double ccart[3], double invert_primitive_basis[3][3], double clattice[3]){
	vecmul(ccart, invert_primitive_basis, clattice);
}

// convert lattice coordinates to cartesian coordinates
void lattice2cartesian(int clattice[3], double primitive_basis[3][3], double ccart[3]){
	for (int dim_idx = 0; dim_idx < 3; dim_idx++){ // vecmul
		ccart[dim_idx] = 
			primitive_basis[dim_idx][0] * (double)clattice[0]
			+ primitive_basis[dim_idx][1] * (double)clattice[1]
			+ primitive_basis[dim_idx][2] * (double)clattice[2];
	}
}

// round floating-point val towards nearest integer in direction of target
int round_towards(double val, int target)
{
	// ENHANCE: do using math.h rounding modes
	if (target >= val)
		return (int) ceil(val);
	else
		return (int) floor(val);
}

// converts lattice basis vectors to unit cell parameters
void primitive_basis2ucell_params(double primitive_basis[3][3], double ucell_params[6]) // primitive_basis = basis vectors (rows/first index), ucell_params = unit cell parameters
{
	double rad2deg = 180.0/PI; // radians to degrees conversion factor
	// a b c - magnitude of basis0 basis1 basis2 vectors
	// ENHANCE: if a vector was primitive_basis[0][*], then this could be done with fdot(u,u) and mag(u)
	// TODO: flip indices of primitive_basis
	ucell_params[0] = sqrt(primitive_basis[0][0]*primitive_basis[0][0] + primitive_basis[1][0]*primitive_basis[1][0] + primitive_basis[2][0]*primitive_basis[2][0]);
	ucell_params[1] = sqrt(primitive_basis[0][1]*primitive_basis[0][1] + primitive_basis[1][1]*primitive_basis[1][1] + primitive_basis[2][1]*primitive_basis[2][1]);
	ucell_params[2] = sqrt(primitive_basis[0][2]*primitive_basis[0][2] + primitive_basis[1][2]*primitive_basis[1][2] + primitive_basis[2][2]*primitive_basis[2][2]);
	// gamma - angle between basis0 and basis1 = arccos(fdot(basis0, basis1) / (mag(basis0) * mag(basis1))); from cos(theta) = fdot(a,b) / (mag(a)*mag(b))
	ucell_params[5] = primitive_basis[0][0]*primitive_basis[0][1] + primitive_basis[1][0]*primitive_basis[1][1] + primitive_basis[2][0]*primitive_basis[2][1];
	ucell_params[5] = ucell_params[5]/(ucell_params[0]*ucell_params[1]);
	ucell_params[5] = rad2deg*acos(ucell_params[5]);
	// beta - angle between basis0 and basis2
	ucell_params[4] = primitive_basis[0][0]*primitive_basis[0][2] + primitive_basis[1][0]*primitive_basis[1][2] + primitive_basis[2][0]*primitive_basis[2][2];
	ucell_params[4] = ucell_params[4]/(ucell_params[0]*ucell_params[2]);
	ucell_params[4] = rad2deg*acos(ucell_params[4]);
	// alpha - angle between basis1 and basis2
	ucell_params[3] = primitive_basis[0][2]*primitive_basis[0][1] + primitive_basis[1][2]*primitive_basis[1][1] + primitive_basis[2][2]*primitive_basis[2][1];
	ucell_params[3] = ucell_params[3]/(ucell_params[1]*ucell_params[2]);
	ucell_params[3] = rad2deg*acos(ucell_params[3]);
	// ENHANCE: double arithmetic leads to imprecise values
	return;
}

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

void adjust_pbc(int* u, int* v, int* w, struct SimulationEnv* se) // should be lattice coordinates
{
	// TODO: allow to turn off pbc in a direction
	// if exceeds boundary, either:
	// treat as occupied site - of what composition? atoms will scale box walls like adatoms
	// *treat as not a site - doesn't contribute to energy, can't be transitioned to
	// 		if not a site, return -2 for adjust_pbc, skip atom_at+findzone, don't collect energy or create transition

	// x y z in lattice coordinates

	if (*u < se->simbox_limits_lat[0][0])
		*u += se->lat_range[0];
	if (*u >= se->simbox_limits_lat[0][1])
		*u -= se->lat_range[0];

	if (*v < se->simbox_limits_lat[1][0])
		*v += se->lat_range[1];
	if (*v >= se->simbox_limits_lat[1][1])
		*v -= se->lat_range[1];

	if (*w < se->simbox_limits_lat[2][0])
		*w += se->lat_range[2];
	if (*w >= se->simbox_limits_lat[2][1])
		*w -= se->lat_range[2];

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
	*zone_u = (int) (((double) se->zone_count_u / (se->lat_range[0])) * (u - se->simbox_limits_lat[0][0]));
	*zone_v = (int) (((double) se->zone_count_v / (se->lat_range[1])) * (v - se->simbox_limits_lat[1][0]));
	*zone_w = (int) (((double) se->zone_count_w / (se->lat_range[2])) * (w - se->simbox_limits_lat[2][0]));

	return;
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
