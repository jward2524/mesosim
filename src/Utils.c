#include "Utils.h" // includes State.h
#include "Vector.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double ssr;

// for simulation box
static double lhs[6];
static int translation_vector[6][3];

// six sides of box
double normal_cart[6][3] = {
    {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}, {-1., 0., 0.}, {0., -1., 0.}, {0., 0., -1.},
};

const double DEFAULT_EPSILON = 1e-3;

static int int_check(double fvalue, int ireference, double epsilon);
static int fact(int n);

// checks if double is within range of an integer - returns 1 for true, 0 for false
static int int_check(double fvalue, int ireference, double epsilon)
{
    return fabs(fvalue - (double)ireference) < epsilon;
}

/**
 * @brief returns a random decimal on the interval [0,1)
 *
 * @return double random value
 */
double drand(void) { return (double)rand() / RAND_MAX; }

/**
 * @brief Get the corresponding index in atom_env
 * 
 * @param offset_idx index in atom.neighbor_arr and se->transition_vectors
 * @param atom_type type of 'main' atom (interchangable with neighbor_type)
 * @param neighbor_type type of neighbor atom wrt 'main'
 * @param se pointer to SimulationEnv state variable
 * @return int index in atom_env
 */
int get_env_index(int offset_idx, int atom_type, int neighbor_type, struct SimulationEnv *se)
{
    // find nearest-neighbor shell
    int nn_level = 0;
    int diff = offset_idx;
    for (int j = 0; j < se->num_nn_levels; j++) {
        diff = diff - se->atoms_per_nn_level[j];
        if (diff < 0) {
            nn_level = j;
            break;
        }
    }

    int bond_idx = get_bond_index(atom_type, neighbor_type, se->num_elements);
    int env_idx = nn_bondidx_2_envidx(nn_level, bond_idx, se->num_bond_types);
    return env_idx;
}

// gets the index in atom_env, nnE arrays (nearest_neighbor - bond_type combo)
int nn_bondidx_2_envidx(int nn, int bond_idx, int num_bond_types)
{
    // nn - nearest neighbor level minus 1 (0 for 1st nn, 1 for 2nd nn, etc.)
    return nn * num_bond_types + bond_idx;
}

int get_bond_index(int a, int b, int num_elements)
{
    // aa, ab, ac; bb, bc; cc [num_elements=3]
    // 00, 01, 02; 11, 12; 22
    // assume 0-indexed
    int first, second;

    // larger number (later element) is second
    if (a < b) {
        first = a;
        second = b;
    } else {
        first = b;
        second = a;
    }

    // a=1, b=2 -> (1*3)+(2-1)=4
    return (first * num_elements) + (second - first);
}

// calculate the number of bond types
int get_num_bond_types(int num_elements)
{
    return fact(num_elements + 2 - 1) / (fact(2) * fact(num_elements - 1));
}

// returns the factorial of n
static int fact(int n)
{
    switch (n) {
    case 0:
        return 1;
    case 1:
        return 1;
    case 2:
        return 2;
    case 3:
        return 6;
    case 4:
        return 24;
    case 5:
        return 120;
    case 6:
        return 720;
    case 7:
        return 5040;
    default:
        fprintf(stderr, "Factorial is too large (max n is 7): %d", n);
        return -1;
    }
}

void lattice2int(double fcoords[3], int coords[3], double epsilon)
{
    for (int dim_idx = 0; dim_idx < 3; dim_idx++) {
        double u = fcoords[dim_idx];
        int comp = (int)round(u);
        int res = int_check(u, comp, epsilon);
        assert(res == 1);
        coords[dim_idx] = comp;
    }
}

// convert a site from cartesian coordinates to lattice coordinates
void cartesian2lattice_site(double ccart[3], double invert_primitive_basis[3][3], int clattice[3])
{
    double fclattice[3];
    vecmul(ccart, invert_primitive_basis, fclattice);
    lattice2int(fclattice, clattice, DEFAULT_EPSILON);
}

void cartesian2lattice(double ccart[3], double invert_primitive_basis[3][3], double clattice[3])
{
    vecmul(ccart, invert_primitive_basis, clattice);
}

// convert lattice coordinates to cartesian coordinates
void lattice2cartesian(int clattice[3], double primitive_basis[3][3], double ccart[3])
{
    for (int dim_idx = 0; dim_idx < 3; dim_idx++) { // vecmul
        ccart[dim_idx] = primitive_basis[dim_idx][0] * (double)clattice[0] +
                         primitive_basis[dim_idx][1] * (double)clattice[1] +
                         primitive_basis[dim_idx][2] * (double)clattice[2];
    }
}

// round floating-point val towards nearest integer in direction of target
int round_towards(double val, int target)
{
    // ENHANCE: do using math.h rounding modes
    if (target >= val)
        return (int)ceil(val);
    else
        return (int)floor(val);
}

// converts lattice basis vectors to unit cell parameters
void primitive_basis2ucell_params(
    double primitive_basis[3][3],
    double ucell_params[6]) // primitive_basis = basis vectors (rows/first index), ucell_params =
                            // unit cell parameters
{
    double rad2deg = 180.0 / PI; // radians to degrees conversion factor
    // a b c - magnitude of basis0 basis1 basis2 vectors
    // ENHANCE: if a vector was primitive_basis[0][*], then this could be done with fdot(u,u) and
    // mag(u)
    // TODO: flip indices of primitive_basis
    ucell_params[0] = sqrt(primitive_basis[0][0] * primitive_basis[0][0] +
                           primitive_basis[1][0] * primitive_basis[1][0] +
                           primitive_basis[2][0] * primitive_basis[2][0]);
    ucell_params[1] = sqrt(primitive_basis[0][1] * primitive_basis[0][1] +
                           primitive_basis[1][1] * primitive_basis[1][1] +
                           primitive_basis[2][1] * primitive_basis[2][1]);
    ucell_params[2] = sqrt(primitive_basis[0][2] * primitive_basis[0][2] +
                           primitive_basis[1][2] * primitive_basis[1][2] +
                           primitive_basis[2][2] * primitive_basis[2][2]);
    // gamma - angle between basis0 and basis1 = arccos(fdot(basis0, basis1) / (mag(basis0) *
    // mag(basis1))); from cos(theta) = fdot(a,b) / (mag(a)*mag(b))
    ucell_params[5] = primitive_basis[0][0] * primitive_basis[0][1] +
                      primitive_basis[1][0] * primitive_basis[1][1] +
                      primitive_basis[2][0] * primitive_basis[2][1];
    ucell_params[5] = ucell_params[5] / (ucell_params[0] * ucell_params[1]);
    ucell_params[5] = rad2deg * acos(ucell_params[5]);
    // beta - angle between basis0 and basis2
    ucell_params[4] = primitive_basis[0][0] * primitive_basis[0][2] +
                      primitive_basis[1][0] * primitive_basis[1][2] +
                      primitive_basis[2][0] * primitive_basis[2][2];
    ucell_params[4] = ucell_params[4] / (ucell_params[0] * ucell_params[2]);
    ucell_params[4] = rad2deg * acos(ucell_params[4]);
    // alpha - angle between basis1 and basis2
    ucell_params[3] = primitive_basis[0][2] * primitive_basis[0][1] +
                      primitive_basis[1][2] * primitive_basis[1][1] +
                      primitive_basis[2][2] * primitive_basis[2][1];
    ucell_params[3] = ucell_params[3] / (ucell_params[1] * ucell_params[2]);
    ucell_params[3] = rad2deg * acos(ucell_params[3]);
    // ENHANCE: double arithmetic leads to imprecise values
    return;
}

// XXX: only used for re-deposition
void get_system_rw_radius(struct SimulationEnv *se)
{
    int ss;

    // find minimum axial distance to system edge

    ss = se->system_size_x;
    if (se->system_size_y < ss)
        ss = se->system_size_y;
    if (se->system_size_z < ss)
        ss = se->system_size_z;

    ssr = (double)ss / 2.;
    ssr = ssr - 5.;

    return;
}

void adjust_pbc(int *u, int *v, int *w, struct SimulationEnv *se) // should be lattice coordinates
{
    // TODO: allow to turn off pbc in a direction
    // if exceeds boundary, either:
    // treat as occupied site - of what composition? atoms will scale box walls like adatoms
    // *treat as not a site - doesn't contribute to energy, can't be transitioned to
    // 		if not a site, return -2 for adjust_pbc, skip atom_at+findzone, don't collect energy
    // or create transition

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

// finds the zone indices xy yz zz that correspond to the lattice coordinates xxx yyy zzz
void findzone(int *zone_u, int *zone_v, int *zone_w, int u, int v, int w, struct SimulationEnv *se)
{
    // *z are pointers to return indices of the zone, *** are lattice coordinates
    // normalize coordinates to the sblimits, then find which zone
    // (zones / extent) * adjusted_coordinate
    *zone_u =
        (int)(((double)se->zone_count_u / (se->lat_range[0])) * (u - se->simbox_limits_lat[0][0]));
    *zone_v =
        (int)(((double)se->zone_count_v / (se->lat_range[1])) * (v - se->simbox_limits_lat[1][0]));
    *zone_w =
        (int)(((double)se->zone_count_w / (se->lat_range[2])) * (w - se->simbox_limits_lat[2][0]));

    return;
}

// check if coordinate passed periodic boundary conditions, and translate if it did
void check_pbc(int *u, int *v, int *w, double basis[3][3])
{
    int coords_lat[3] = {*u, *v, *w};
    // double fcoords_lat[3] = {(double) *x, (double) *y, (double) *z};
    int rhs;
    double coords_cart[3];
    lattice2cartesian(coords_lat, basis, coords_cart);
    for (int side = 0; side < 6; side++) {
        // rhs contains the point being tested (xyz)
        // negative normal gives negative rhs
        rhs = (int)fdot(normal_cart[side], coords_cart);

        // example conditions to be in region
        // positive normal: 0 (plane) < 2 (coord)
        // negative normal: -128 (plane) < -127 (coord)
        // non-origin sides (3-5) should reroute to 0
        bool cond;
        if (side >= 3)
            cond = lhs[side] > rhs;
        else
            cond = lhs[side] >= rhs;

        if (cond) {
            pbc_translate(coords_lat, translation_vector[side]);
            *u = coords_lat[0];
            *v = coords_lat[1];
            *w = coords_lat[2];
        }
    }
}

// translate a point according to boundary conditions
void pbc_translate(int coords_lat[3], int translation_vector_l[3])
{
    ivecsum(coords_lat, translation_vector_l, coords_lat);
}

/**
 * @brief Get atom type from a name
 * 
 * @param atom_name 
 * @param atom_names array of atom names, from input file usually
 * @param atom_names_cnt length of atom_names array
 * @param atom_type output variable pointer
 * @return int 
 */
int get_type_from_name(char *atom_name, char **atom_names, int atom_names_cnt, unsigned char *atom_type)
{
    int found = 0;
    for (int i = 0; i < atom_names_cnt; i++) {
        int str_match = strncmp(atom_name, atom_names[i], strlen(atom_names[i])) == 0;
        if (str_match) {
            *atom_type = (unsigned char) i;
            found = 1;
            return 0;
        }
    }
    if (found == 0) {
        fprintf(stderr, "Atom name %s not found in list of atom names\n", atom_name);
        return 1;
    } else {
        fprintf(stderr, "Atom name found but didn't return earlier, something went wrong");
        return 1;
    }
}
