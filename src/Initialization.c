#include "Initialization.h"
#include "Atoms.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Utils.h"
#include "Vector.h"
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double normal_x, normal_y, normal_z;

// for simulation box
static double lhs[6];
static int translation_vector[6][3];

// direction of possible atom jumps for each crystal lattice type
const LatticeVector BCC_OFFSET[8] = { // not normalized, in lattice coordinates
    {-1, -1, -1}, {0, 0, -1}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {1, 1, 1}};

const LatticeVector FCC_OFFSET[12] = {{0, 1, -1}, {1, 0, -1}, {0, 0, -1}, {1, 0, 0},
                                      {-1, 0, 0}, {0, 1, 0},  {0, -1, 0}, {0, 0, 1},
                                      {-1, 1, 0}, {1, -1, 0}, {-1, 0, 1}, {0, -1, 1}};

const LatticeVector SC_OFFSET[6] = {{0, 0, -1}, {0, 0, 1}, {0, 1, 0},
                                    {0, -1, 0}, {1, 0, 0}, {-1, 0, 0}};

const LatticeVector SC_OFFSET_2[12] = {{1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},
                                       {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
                                       {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}};

const LatticeVector FCC_OFFSET_2[6] = {{1, -1, 1},  {-1, 1, 1}, {-1, -1, 1},
                                       {-1, 1, -1}, {1, 1, -1}, {1, -1, -1}};

const LatticeVector BCC_OFFSET_2[6] = {{1, 0, 1},   {-1, 0, -1}, {1, 1, 0},
                                       {-1, -1, 0}, {0, 1, 1},   {0, -1, -1}};

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

    cross(a, b, n);

    normal_x = n[0];
    normal_y = n[1];
    normal_z = n[2];

    nmag = sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z);

    normal_x = normal_x / nmag;
    normal_y = normal_y / nmag;
    normal_z = normal_z / nmag;

    if (normal_z < 0) {
        normal_x *= -1.0;
        normal_y *= -1.0;
        normal_z *= -1.0;
    }

    return;
}

// updates [iv, iy; primitive_basis, ucell_params, Atoms' cart_coords?; zone_arr, rmat; normal_x,
// normal_y, normal_z; num_transition_vectors, se->transition_vectors, se->opposite_tvectors; zi*,
// zi*shift, *sh], rate_cnt, transition_cnt, atom_cnt, frequency_sum, elapsed_stime, overpotential,
// next_log_checkpoint
void initialize_simulation_variables(struct SimulationState *ss, struct SimulationEnv *se)
{
    // holy shit, when not =DISSOLUTION[=1], it breaks some shit heavy
    if (!se->dissolution) {
        se->is_soluble = (bool *)malloc((size_t)se->num_elements * sizeof(bool));
        for (int i = 0; i < se->num_elements; i++) {
            se->is_soluble[i] = false;
        }
    }

    // finish initializing structs

    // max_rates = combination of contributors and number of nn bond types
    // (and evaporation and final configuration count)
    // each contributor can be any kind of bond:
    // (1st nn bonds) ^ (1st nn contributors) + (2st nn bonds) ^ (2st nn contributors)
    // pow(base, exponent)
    se->max_rates = 0;
    for (int i = 0; i < se->num_nn_levels; i++) {
        se->max_rates += (long int)pow(se->num_bond_types, se->atoms_per_nn_level[i]);
    }
    se->max_rates = (long int)((1 + se->dissolution) * se->max_rates);

    se->max_atoms = se->lat_range[0] * se->lat_range[1] * se->lat_range[2];
    se->max_transitions = ((se->num_transition_vectors + se->dissolution) * se->max_atoms) + 10;

    ss->atom_arr = (Atom **)malloc((size_t)se->max_atoms * sizeof(Atom *));
    ss->rate_arr = (Rate *)malloc((size_t)se->max_rates * sizeof(Rate));
    ss->transition_arr = (Transition **)malloc((size_t)se->max_transitions * sizeof(Transition *));

    // null pointer checks
    if (!ss->atom_arr) {
        perror("Couldn't allocate memory for atom array");
        clean_and_error(errno);
    }
    if (!ss->rate_arr) {
        perror("Couldn't allocate memory for rate array");
        clean_and_error(errno);
    }
    if (!ss->transition_arr) {
        perror("Couldn't allocate memory for transition array");
        clean_and_error(errno);
    }

    ss->transition_probability.rate_arr_index =
        (long *)malloc((size_t)se->max_rates * sizeof(long));
    ss->transition_probability.lbound = (double *)malloc((size_t)se->max_rates * sizeof(double));
    ss->transition_probability.ubound = (double *)malloc((size_t)se->max_rates * sizeof(double));

    if (!ss->transition_probability.rate_arr_index) {
        perror("Couldn't allocate memory for transition probability rate_arr_index");
        clean_and_error(errno);
    }
    if (!ss->transition_probability.lbound) {
        perror("Couldn't allocate memory for transition probability lbound");
        clean_and_error(errno);
    }
    if (!ss->transition_probability.ubound) {
        perror("Couldn't allocate memory for transition probability ubound");
        clean_and_error(errno);
    }

    ss->rate_cnt = 0; // initialize global transition variables
    ss->transition_cnt = 0;
    ss->atom_cnt = 0; // initialize global atom variables
    // XXX: commented code, never used
    // current_iteration = 0; //not needed if only running 1 simulation at a time
    ss->frequency_sum = 0.0;

    // TODO: allow for reading simulation variables from intermediate xyz file
    ss->elapsed_stime = 0.0;

    ss->total_atoms_dissolved = 0;

    // iteration count
    ss->iter = 0;

    ss->overpotential = se->initial_overpotential;

    srand(se->rand_seed);

    se->centroid[0] = 0.;
    se->centroid[1] = 0.;
    se->centroid[2] = 0.;

    // // XXX: redundant
    // // next_log_checkpoint is initialized to one log_interval_step step
    // if (ls->analysis_type == REGULAR_TIME_INTERVALS) // TODO: reconsider what is happening here
    // 	ls->next_log_checkpoint = ls->next_log_checkpoint;	//do we want this to be true? -
    // overwrites what was in the input file else if (ls->analysis_type == LN_TIME_INTERVALS)
    // ls->next_log_checkpoint = ls->next_log_checkpoint;

    // first, remove any atoms that may exist
    // [ ]: why would atom_cnt not be zero????
    while (ss->atom_cnt != 0) // TODO: start here
        kill_atom(ss->atom_cnt - 1, ss, se);

    return;
}

void initialize_initial_structure(struct SimulationState *ss, struct SimulationEnv *se,
                                  struct LoggingState *ls)
{
    // index represents geometry, from macros
    switch (se->geometry) {
    case GEOMETRY_FLAT_SHEET: // flat plane
        initialize_flat_sheet(ss, se);
        break;

    case GEOMETRY_CLUSTER:
        initialize_spherical_cluster(ss, se);
        break;

    case GEOMETRY_FROM_FILE:
        initialize_from_file(ss, se, ls);
        break;
    }

    // optimizes the atoms added in the initialization routines??
    // check_system(ss, se); // XXX
    // organize(ss->atom_arr, ss->atom_cnt, se->primitive_basis);

    return;
}

/********************************************************************************/
/********************************************************************************/
// zi* are the number of zones in that dimension, zi*shift is for bit shifting to find which zone a
// lattice coordinate corresponds to?
void get_shifts(struct SimulationEnv *se)
{ // updates zi*, zi*shift, *sh
    int temp1;
    se->zone_count_u = TTS;
    se->zone_count_v = TTS;
    se->zone_count_w = TTS;

    temp1 = se->zone_count_u;
    se->zixshift = 0;
    while (temp1 > 1) {
        ++se->zixshift;
        temp1 = temp1 / 2;
    }

    temp1 = se->zone_count_v;
    se->ziyshift = 0;
    while (temp1 > 1) {
        ++se->ziyshift;
        temp1 = temp1 / 2;
    }

    temp1 = se->zone_count_w;
    se->zizshift = 0;
    while (temp1 > 1) {
        ++se->zizshift;
        temp1 = temp1 / 2;
    }

    temp1 = se->system_size_x;
    se->ssxshift = 0;
    while (temp1 > 1) {
        ++se->ssxshift;
        temp1 = temp1 / 2;
    }

    temp1 = se->system_size_y;
    se->ssyshift = 0;
    while (temp1 > 1) {
        ++se->ssyshift;
        temp1 = temp1 / 2;
    }

    temp1 = se->system_size_z;
    se->sszshift = 0;
    while (temp1 > 1) {
        ++se->sszshift;
        temp1 = temp1 / 2;
    }
    // never used, just left and right shift with se->zixshift and se->ssxshift in findzone()
    se->xsh = se->zixshift - se->ssxshift;
    se->ysh = se->ziyshift - se->ssyshift;
    se->zsh = se->zizshift - se->sszshift;
    // TODO: more shifts in zones than in system? se->zone_count_u > se->system_size_x
    return;
}

// updates rmat
// supposedly for viewing
void set_default_orientation(int lattice_type, double rmat[3][3])
{
    double axis[3], pnormal[3], axis_mag;
    double zero_point[3] = {0., 0., 0.}, a_a[3];
    double zaxis[3] = {0., 0., 1.};
    double zeroa[3] = {0., 0., 0.};
    double spin_ax[3];
    double vec_angle;

    // organize(atom_arr, atom_cnt, basis); // atom_cnt ='d 0

    switch (lattice_type) {
    case FCC:
        axis[0] = -1.0; // {1,1,1}
        axis[1] = 1.0;
        axis[2] = 1.0;
        break;

    case BCC:
        axis[0] = 0.0; // {1,1,1}
        axis[1] = 1.0;
        axis[2] = -1.0;
        break;

    case SC:
    default:
        axis[0] = 0.0; // {1,1,1}
        axis[1] = 0.0;
        axis[2] = 1.0;
        break;
    }

    axis_mag = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);

    // direction cosines of axis normal to plane

    pnormal[0] = axis[0] / axis_mag;
    pnormal[1] = axis[1] / axis_mag;
    pnormal[2] = axis[2] / axis_mag;

    if ((pnormal[0] == zaxis[0]) && (pnormal[1] == zaxis[1]) && (pnormal[1] == zaxis[1])) {
        identity2(rmat);
        // removed rotation_notify_flag
        return;
    }

    // orient
    vecdif(pnormal, zero_point, a_a);
    if (((a_a[0] == 0.) && (a_a[1] == 0.)) || (magnitude(a_a) == 0.))
        return;

    normto(a_a, zaxis, spin_ax);
    vec_angle = -0.0174533 * vangle(zaxis, zeroa, a_a);

    rotmata(spin_ax, vec_angle, rmat);
    transpose(rmat);

    // rotation_notify_flag removed

    return;
}

// updates zone (array) based on zi* (zone sizes?), initializes offset to -1
void initialize_zones(Zone zone_arr[ZONES_IN_X][ZONES_IN_Y][ZONES_IN_Z], struct SimulationEnv *se)
{
    int i, j, k;

    for (i = 0; i < se->zone_count_u; ++i)
        for (j = 0; j < se->zone_count_v; ++j)
            for (k = 0; k < se->zone_count_w; ++k)
                zone_arr[i][j][k].offset = -1;
    return;
}

// initializes primitve_basis, ucell_params, ss*
void initialize_lattice_geometry(struct SimulationEnv *se)
{
    // Initializes the generic lattice geometry to be simple cubic (i.e., a=1, b=1, c=1, alpha = 90,
    // beta = 90, gamma = 90)

    // int i,j;

    // for (i=0;i<3;++i)
    // 	for (j=0;j<3;++j)
    // 		if (i == j) primitive_basis[i][j] = 1.0; else primitive_basis[i][j] = 0.0;

    double pb[3][3] = {
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.},
    };
    memcpy(se->primitive_basis, pb, 3 * 3 * sizeof(double));

    inver(se->primitive_basis, se->invert_primitive_basis);
    primitive_basis2ucell_params(se->primitive_basis, se->ucell_params);

    se->system_size_x = 1;
    se->system_size_y = 1;
    se->system_size_z = 1;

    return;
}

/********************************************************************************/
/********************************************************************************/
// updates num_transition_vectors, [se->transition_vectors, se->opposite_tvectors]
void initialize_neighbor_offsets(struct SimulationEnv *se)
{
    int i;
    int fcc_opposite[12] = {11, 10, 7, 4, 3, 6, 5, 2, 9, 8, 1, 0};
    int sc_opposite[6] = {1, 0, 3, 2, 5, 4};
    int bcc_opposite[8] = {7, 6, 3, 2, 5, 4, 1, 0};

    const LatticeVector *vectors;
    const LatticeVector *vectors2;
    const int *opposite_vectors;
    int extra;

    // assumes only first-nearest neighbor transitions
    // but can energy contributors can be second-nearest neighbors
    switch (se->lattice_type) {
    case FCC:
        se->num_transition_vectors = 12;
        extra = 6;
        vectors = FCC_OFFSET;
        vectors2 = FCC_OFFSET_2;
        opposite_vectors = fcc_opposite;
        break;

    case SC:
        se->num_transition_vectors = 6;
        extra = 12;
        vectors = SC_OFFSET;
        vectors2 = SC_OFFSET_2;
        opposite_vectors = sc_opposite;
        break;

    case BCC:
        se->num_transition_vectors = 8;
        extra = 6;
        vectors = BCC_OFFSET;
        vectors2 = BCC_OFFSET_2;
        opposite_vectors = bcc_opposite;
        break;
    }

    int transition_length;
    se->atoms_per_nn_level = (int *)malloc((size_t)se->num_nn_levels * sizeof(int));
    se->atoms_per_nn_level[0] = se->num_transition_vectors;
    if (se->num_nn_levels == 1) {
        se->num_energy_contributors = se->num_transition_vectors;
    } else if (se->num_nn_levels == 2) {
        se->atoms_per_nn_level[1] = extra;
        se->num_energy_contributors = se->num_transition_vectors + extra;
    }
    transition_length = se->num_energy_contributors;

    se->transition_vectors =
        (LatticeVector *)malloc((size_t)transition_length * sizeof(LatticeVector));
    if (!se->transition_vectors) {
        perror("Couldn't allocate memory for jump offset array");
        clean_and_error(errno);
    }

    // opposite_tvectors always only 1st nn
    se->opposite_tvectors =
        (int *)malloc((size_t)se->num_transition_vectors * sizeof(*se->opposite_tvectors));
    if (!se->opposite_tvectors) {
        perror("Couldn't allocate memory for jump offset array");
        clean_and_error(errno);
    }

    for (i = 0; i < se->num_transition_vectors; ++i) {
        se->transition_vectors[i].dx = vectors[i].dx;
        se->transition_vectors[i].dy = vectors[i].dy;
        se->transition_vectors[i].dz = vectors[i].dz;
        se->opposite_tvectors[i] = opposite_vectors[i];
    }

    if (se->num_nn_levels == 2) {
        for (i = 0; i < extra; i++) {
            se->transition_vectors[se->num_transition_vectors + i].dx = vectors2[i].dx;
            se->transition_vectors[se->num_transition_vectors + i].dy = vectors2[i].dy;
            se->transition_vectors[se->num_transition_vectors + i].dz = vectors2[i].dz;
        }
    }

    return;
}

/********************************************************************************/
/********************************************************************************/
// updates primitive_basis, ucell_params, Atoms' cart_coords? to match lattice type
void set_primitive_basis(struct SimulationEnv *se) // lattice_type = crystal structure type
{
    switch (se->lattice_type) {
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

void initialize_flat_sheet(struct SimulationState *ss, struct SimulationEnv *se)
{
    double rand;

    int mid_w = se->simbox_limits_lat[2][0] +
                ((se->simbox_limits_lat[2][1] - se->simbox_limits_lat[2][0]) / 2);
    int half_thickness = se->sheet_thickness / 2;

    for (int k = mid_w - half_thickness; k < mid_w + half_thickness; ++k) {
        for (int i = se->simbox_limits_lat[0][0]; i < se->simbox_limits_lat[0][1]; ++i) {
            // loop through x and y (two lattice directions)
            for (int j = se->simbox_limits_lat[1][0]; j < se->simbox_limits_lat[1][1]; ++j) {
                rand = drand();

                // TODO: make into fxn
                double bar = 0;
                int type = -1;
                do {
                    type++;
                    bar += se->substrate_composition[type];
                } while (bar < rand);
                add_atom(i, j, k, (unsigned char)type, NORMAL, ss, se);
            }
        }
    }
    return;
}

/********************************************************************************/
/********************************************************************************/
// ENHANCE: currently adds one extra atom to radius - remove it
// radius of cluster in number of atoms (nearest-neighbor distances)
void initialize_spherical_cluster(struct SimulationState *ss, struct SimulationEnv *se)
{
    double center_cart[3];   // cartesian/orthogonal coordinates of center point
    int center_lattice[3];   // lattice coordinates of center point
    int atom_pos_lattice[3]; // atom position in lattice coords
    double atom_pos_cart[3];

    double random_num; // random number

    double radius_cart; // radius of cluster, in cartesian units

    double dist;

    // for when ss* represents a prism cell
    // center of the cluster is the halfway point
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
    for (int dim = 0; dim < 3; dim++) {
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
    double bbcorners_cart[8][3] = {{bblimits_cart[0][0], bblimits_cart[1][0], bblimits_cart[2][0]},
                                   {bblimits_cart[0][0], bblimits_cart[1][0], bblimits_cart[2][1]},
                                   {bblimits_cart[0][0], bblimits_cart[1][1], bblimits_cart[2][0]},
                                   {bblimits_cart[0][0], bblimits_cart[1][1], bblimits_cart[2][1]},
                                   {bblimits_cart[0][1], bblimits_cart[1][0], bblimits_cart[2][0]},
                                   {bblimits_cart[0][1], bblimits_cart[1][0], bblimits_cart[2][1]},
                                   {bblimits_cart[0][1], bblimits_cart[1][1], bblimits_cart[2][0]},
                                   {bblimits_cart[0][1], bblimits_cart[1][1], bblimits_cart[2][1]}};

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
            clean_and_error(EXIT_FAILURE);
        }
        if (bblimits_cart[dim_idx][1] > (int)(2 * center_cart[dim_idx])) {
            printf("ERROR! Spherical cluster passes through periodic boundary conditions\n");
            clean_and_error(EXIT_FAILURE);
        }
    }

    // iterates through the bounding box of the sphere to identify positions in cluster
    // ENHANCE: only 52% of loops will be successful - make it more efficient
    for (int u = bblimits_lattice[0][0]; u <= bblimits_lattice[0][1]; u++) {
        for (int v = bblimits_lattice[1][0]; v <= bblimits_lattice[1][1]; v++) {
            for (int w = bblimits_lattice[2][0]; w <= bblimits_lattice[2][1]; w++) {
                // convert i, j, k to cartesian coordinates
                atom_pos_lattice[0] = u;
                atom_pos_lattice[1] = v;
                atom_pos_lattice[2] = w;

                // to cartesian coordinates
                lattice2cartesian(atom_pos_lattice, se->primitive_basis, atom_pos_cart);
                // distance to center
                dist = (atom_pos_cart[0] - center_cart[0]) * (atom_pos_cart[0] - center_cart[0]) +
                       (atom_pos_cart[1] - center_cart[1]) * (atom_pos_cart[1] - center_cart[1]) +
                       (atom_pos_cart[2] - center_cart[2]) * (atom_pos_cart[2] - center_cart[2]);

                if (dist <= (radius_cart * radius_cart)) {
                    // particle is in bounds
                    random_num = drand();
                    // determining composition of atom to be placed
                    double bar = 0;
                    int type = -1;
                    do {
                        type++;
                        bar += se->substrate_composition[type];
                    } while (bar < random_num);
                    add_atom(u, v, w, (unsigned char)type, NORMAL, ss, se);
                }
            }
        }
    }
    return;
}

/********************************************************************************/
/********************************************************************************/
// get atom positions from a file
void initialize_from_file(struct SimulationState *ss, struct SimulationEnv *se,
                          struct LoggingState *ls)
{
    char *file_ext = strrchr(se->atoms_filename, '.') + 1;
    if (strcmp(file_ext, "xyz") != 0) {
        fprintf(stderr, "File extension %s not recognized as atom input file", file_ext);
        clean_and_error(EXIT_FAILURE);
    }

    FILE *atom_file = fopen(se->atoms_filename, "r");
    if (!atom_file) {
        printf("ERROR! Couldn't open output file %s\n", se->atoms_filename);
        fprintf(stderr, "Couldn't open file %s: %s\n", se->atoms_filename, strerror(errno));
        clean_and_error(errno);
    }

    process_xyz_file(ls->sim_log, atom_file, ss, se, ls);
    // TODO: compatability check between parameters in comment line of xyz file and input file
    return;
}

// int simbox_limits_lat[3][2]; // lattice limits of simulation box in each dimension - for zones
// double system_size_x, double system_size_y, double system_size_z)
void initialize_simulation_box(struct SimulationEnv *se)
{
    // assuming simulation box/prism
    // system size in cartesian units [nearest-neighbor (or some other lattice-based) units in
    // cartesian grid] system size of 128 -> 0 to 127, 128th->0 need to define box in terms of
    // lattice vectors 6 planes, of form dot(normal, point on plane) = dot(normal, [x,y,z of point
    // to test]) normal=(1,0,0); point on plane=(128,0,0) -> 128 = x normal of family (1,0,0)
    // (-1,0,0) point on plane of family  (se->system_size_x,0,0) (0,0,0) x=0 y=0 z=0
    // x=se->system_size_x y=se->system_size_y z=se->system_size_z to define a region: x>=0 y>=0
    // z>=0 x<se->system_size_x y<se->system_size_y z<se->system_size_z normals point towards inside
    // of region (keeps the inequality the same)

    // convert into lattice vector form
    // dot(cart2lattice(normal), cart2lattice(point on plane)) {fixed} = dot(cart2lattice(normal),
    // [point to test]) [6 planes x 1] vector = [6 normals x 3 coords] @ [3 x 1] = [6 x 1] vector
    // ([6] vector * [6] sign for comparison + [0 or -1] for counting on the boundary) >= 0

    // each plane also has associated translation vector, cart2lattice(system size * normal)
    // if outside region, translate by translation vector

    double point_cart[6][3] = {
        {0, 0, 0},                 // x
        {0, 0, 0},                 // y
        {0, 0, 0},                 // z
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

    for (int i = 0; i < 3; i++) {
        se->simbox_limits_lat[i][0] = center_lattice[0];
        se->simbox_limits_lat[i][1] = center_lattice[0];
    }
    corners2limits(sbcorners_cart, se->simbox_limits_lat, se->invert_primitive_basis);

    for (int i = 0; i < 3; i++) {
        se->lat_range[i] = se->simbox_limits_lat[i][1] - se->simbox_limits_lat[i][0];
        int vector_lat[3] = {0, 0, 0};
        vector_lat[i] = se->lat_range[i];
        lattice2cartesian(vector_lat, se->primitive_basis, se->simbox_vectors_cart[i]);
    }

    // origin is the minimum in all directions
    int origin_lat[] = {se->simbox_limits_lat[0][0], se->simbox_limits_lat[1][0],
                        se->simbox_limits_lat[2][0]};

    lattice2cartesian(origin_lat, se->primitive_basis, se->simbox_origin_cart);

    // double point_lat[6][3];
    double translation_dist;
    double scaled_normal_cart[3];
    double normal_lat[6][3]; // TODO: when reimplemented, if needs to have larger scope
    for (int side = 0; side < 6; side++) {
        // doesn't contribute; only for sake of seeing normal in lat coordinates
        vecmul(normal_cart[side], se->invert_primitive_basis, normal_lat[side]);

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
        cartesian2lattice_site(scaled_normal_cart, se->invert_primitive_basis,
                               translation_vector[side]);
    }
}

// TODO:
void smallest_translation(double normal_cart_t[3], double inv_basis[3][3])
{
    // first, get base conversion
    double normal_lat[3];
    vecmul(normal_cart_t, inv_basis, normal_lat);

    // if conversion has decimals, invert them and try to convert to integers
    // if conversion fails, remove whole number from inversion and try to invert again
    // repeat, multiplying whole numbers together as you go
    // if all integers are the same, multiply normal_lat by integer
    // if not all the same integers, find Least Common Multiple and multipy normal_lat by LCM
}

// converts cartesian corners to lattice limits along each dimension
// values are rounded towards initial values of limits_lat
void corners2limits(double corners_cart[8][3], int limits_lat[3][2], double inv_basis[3][3])
{
    // convert corners from cartesian coords into atom/lattice coords
    // and find limits in each dimension
    double bbcorners_lattice[8][3];
    // TODO: make data types of cart and lattice coordinates make more sense, and thus the
    // conversion functions find the loop limits from min and max coordinates of corners
    for (int corner_idx = 0; corner_idx < 8; corner_idx++) {
        cartesian2lattice(corners_cart[corner_idx], inv_basis, bbcorners_lattice[corner_idx]);

        for (int dim_idx = 0; dim_idx < 3; dim_idx++) {
            // if coordinate is smaller than lower limit, change limit (round towards center)
            if (bbcorners_lattice[corner_idx][dim_idx] < (double)limits_lat[dim_idx][0])
                limits_lat[dim_idx][0] =
                    round_towards(bbcorners_lattice[corner_idx][dim_idx], limits_lat[dim_idx][0]);

            // if coordinate is largger than upper limit, change limit (round towards center)
            else if (bbcorners_lattice[corner_idx][dim_idx] > (double)limits_lat[dim_idx][1])
                limits_lat[dim_idx][1] =
                    round_towards(bbcorners_lattice[corner_idx][dim_idx], limits_lat[dim_idx][1]);
        }
    }
}

void initialize_simulation(struct SimulationState *ss, struct SimulationEnv *se,
                           struct LoggingState *ls)
{

    // bit shifts for periodic boundary conditions
    // requires: zone_count_uvw, system_size_xyz
    // updates: zixyzshift, ssxyzshift, xyzsh
    get_shifts(se);

    // initialize zones - help figure out which atoms are next to which other atoms
    // requires: zone_count_uvw
    // updates: zone_arr
    initialize_zones(ss->zone_arr, se);

    // requires: lattice_type
    // updates: primitive_basis, invert_primitive_basis, ucell_params
    set_primitive_basis(se);

    // supposedly was only for visualization
    // set_default_orientation(sim_state->atom_arr, sim_state->atom_cnt, sim_env->lattice_type,
    //                         sim_env->rmat, sim_env->primitive_basis);
    // maybe only for visualization
    // get_system_normal(sim_env->primitive_basis);

    // requires: invert_primitive_basis, system_size_xyz
    // updates: simbox_limits_lat, lat_range, simbox_vectors_cart, simbox_origin_cart
    initialize_simulation_box(se);

    /*
    sets jump offsets for given crystal type
    requires: lattice_type, num_nn_levels
    updates: num_transition_vectors, atoms_per_nn_level, num_energy_contributors,
    transition_vectors, opposite_tvectors
    */
    initialize_neighbor_offsets(se);

    // mallocs and sets to zero (or something else) simulation variables
    // requires: dissolution, num_elements, num_bond_types, atoms_per_nn_level, max_atoms,
    // initial_overpotential, lat_range
    // updates: a lot
    initialize_simulation_variables(ss, se);

    // requires: geometry, system_size_xyz, primitive_basis, substrate_composition,
    // simbox_limits_lat, sheet_thickness, cluster_radius
    // updates: atom_arr, atom_cnt, ...
    initialize_initial_structure(ss, se, ls);

    // check_system(sim_state, sim_env); // XXX?

    // default flavor is KMC, if unset
    if (se->flavor == 0)
        se->flavor = FLAVOR_KMC;

    input_logging(ss, se, ls);
}
