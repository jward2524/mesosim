#include "ErrorM.h"
#include "Initialization.h"
#include "TUtils.h"
#include "unity.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>

struct SimulationState *ss;
struct SimulationEnv *se;
struct LoggingState *ls;
FILE *temp_log;

void setUp(void)
{
    initialize_states(&ss, &se, &ls);
    init_temp(&temp_log);
}

void tearDown(void)
{
    // fclose needs to be here in case a test fails
    // set to null to prevent double free from fclose + clean_and_error
    ls->sim_log = NULL;
    fclose(temp_log);
    clean_and_error(0);
}

void test_get_shifts_sets_expected_shifts(void)
{
    // set distinct powers-of-two for zones and system sizes
    se->zone_count_u = 8;
    se->zone_count_v = 4;
    se->zone_count_w = 2;

    se->system_size_x = 8;
    se->system_size_y = 4;
    se->system_size_z = 2;

    get_shifts(se);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->zixshift, "zixshift should be log2(zone_count_u)");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->ziyshift, "ziyshift should be log2(zone_count_v)");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->zizshift, "zizshift should be log2(zone_count_w)");

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->ssxshift, "ssxshift should be log2(system_size_x)");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->ssyshift, "ssyshift should be log2(system_size_y)");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->sszshift, "sszshift should be log2(system_size_z)");

    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->xsh, "xsh should be zixshift - ssxshift (0)");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->ysh, "ysh should be ziyshift - ssyshift (0)");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->zsh, "zsh should be zizshift - sszshift (0)");
}

void test_set_primitive_basis_fcc_sets_expected_vectors(void)
{
    set_primitive_basis(FCC, se);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCXV1, se->primitive_basis[0][0], "fcc basis[0][0]");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCXV2, se->primitive_basis[0][1], "fcc basis[0][1]");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCXV3, se->primitive_basis[0][2], "fcc basis[0][2]");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCYV1, se->primitive_basis[1][0], "fcc basis[1][0]");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCYV2, se->primitive_basis[1][1], "fcc basis[1][1]");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCYV3, se->primitive_basis[1][2], "fcc basis[1][2]");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCZV1, se->primitive_basis[2][0], "fcc basis[2][0]");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCZV2, se->primitive_basis[2][1], "fcc basis[2][1]");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(FCCZV3, se->primitive_basis[2][2], "fcc basis[2][2]");
}

void test_initialize_neighbor_offsets_fcc_sets_expected_values(void)
{
    struct UserInputs inputs = {0};
    inputs.lattice_type = FCC;
    inputs.num_nn_levels = 2;

    initialize_neighbor_offsets(&inputs, se);

    TEST_ASSERT_EQUAL_INT_MESSAGE(12, se->num_transition_vectors, "FCC has 12 first-nn vectors");
    TEST_ASSERT_EQUAL_INT_MESSAGE(18, se->num_energy_contributors,
                                  "FCC with 2 nn levels -> 12+6=18 contributors");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12, se->atoms_per_nn_level[0],
                                  "atoms_per_nn_level[0] should be 12");
    TEST_ASSERT_EQUAL_INT_MESSAGE(6, se->atoms_per_nn_level[1],
                                  "atoms_per_nn_level[1] should be 6");
}

void test_corners2limits_updates_limits_correctly(void)
{
    double corners_cart[8][3] = {{0, 0, 0}, {4, 0, 0}, {0, 4, 0}, {0, 0, 4},
                                 {4, 4, 0}, {4, 0, 4}, {0, 4, 4}, {4, 4, 4}};
    int limits_lat[3][2] = {{2, 2}, {2, 2}, {2, 2}};
    double inv_basis[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    corners2limits(corners_cart, limits_lat, inv_basis);

    TEST_ASSERT_EQUAL_INT_MESSAGE(0, limits_lat[0][0], "x lower limit should become 0");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, limits_lat[0][1], "x upper limit should become 4");

    TEST_ASSERT_EQUAL_INT_MESSAGE(0, limits_lat[1][0], "y lower limit should become 0");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, limits_lat[1][1], "y upper limit should become 4");

    TEST_ASSERT_EQUAL_INT_MESSAGE(0, limits_lat[2][0], "z lower limit should become 0");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, limits_lat[2][1], "z upper limit should become 4");
}

void test_initialize_state_from_input_copies_values(void)
{
    struct UserInputs in = {0};
    in.temperature = 400.0;
    in.initial_overpotential = 0.5;
    in.sim_end_type = SIM_END_BY_STIME;
    in.run_stime = 12.5;
    in.final_iteration = 42;

    initialize_state_from_input(&in, ss);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(400.0, ss->temperature, "temperature copied");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.5, ss->overpotential, "overpotential copied");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_STIME, ss->sim_end_type, "sim_end_type copied");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(12.5, ss->run_stime, "run_stime copied");
    TEST_ASSERT_EQUAL_UINT64_MESSAGE(42, ss->final_iteration, "final_iteration copied");
}

void test_initialize_zones_allocates_and_sets_offsets(void)
{
    // small zone counts to keep allocations tiny
    se->zone_count_u = 2;
    se->zone_count_v = 2;
    se->zone_count_w = 2;

    initialize_zones(&ss->zone_arr, se);

    for (size_t i = 0; i < se->zone_count_u; ++i) {
        for (size_t j = 0; j < se->zone_count_v; ++j) {
            for (size_t k = 0; k < se->zone_count_w; ++k) {
                TEST_ASSERT_EQUAL_INT_MESSAGE(-1, ss->zone_arr[i][j][k].offset,
                                              "zone offset initialized to -1");
            }
        }
    }
}

void test_initialize_simulation_box_and_lat_range(void)
{
    struct UserInputs in = {0};
    in.system_size_x = 8;
    in.system_size_y = 8;
    in.system_size_z = 8;

    // primitive basis must be set before calling initialize_simulation_box
    set_primitive_basis(SC, se);

    initialize_simulation_box(&in, se);

    TEST_ASSERT_EQUAL_INT_MESSAGE(8, se->system_size_x, "system_size_x set");
    TEST_ASSERT_TRUE_MESSAGE(se->lat_range[0] > 0, "lat_range x positive");
    TEST_ASSERT_TRUE_MESSAGE(se->lat_range[1] > 0, "lat_range y positive");
    TEST_ASSERT_TRUE_MESSAGE(se->lat_range[2] > 0, "lat_range z positive");
}

void test_set_default_orientation_sc_returns_identity(void)
{
    double rmat[3][3];
    set_default_orientation(SC, rmat);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.0, rmat[0][0], "rmat[0][0] == 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.0, rmat[0][1], "rmat[0][1] == 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.0, rmat[0][2], "rmat[0][2] == 0");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.0, rmat[1][0], "rmat[1][0] == 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.0, rmat[1][1], "rmat[1][1] == 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.0, rmat[1][2], "rmat[1][2] == 0");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.0, rmat[2][0], "rmat[2][0] == 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.0, rmat[2][1], "rmat[2][1] == 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.0, rmat[2][2], "rmat[2][2] == 1");
}

void test_allocate_simulation_arrays_allocates_expected_buffers(void)
{
    // keep sizes small
    se->zone_count_u = 2;
    se->zone_count_v = 2;
    se->zone_count_w = 2;

    se->max_atoms = 4;
    se->max_rates = 2;
    se->max_transitions = 8;

    allocate_simulation_arrays(ss, se);

    TEST_ASSERT_NOT_NULL_MESSAGE(ss->atom_arr, "atom_arr allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->rate_arr, "rate_arr allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->transition_arr, "transition_arr allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->transition_probability.rate_arr_index,
                                 "transition_probability.rate_arr_index allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->transition_probability.lbound,
                                 "transition_probability.lbound allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->transition_probability.ubound,
                                 "transition_probability.ubound allocated");
}

void test_initialize_env_from_input_sets_seed_and_overpotentials(void)
{
    struct UserInputs in = {0};
    in.system_size_x = 8;
    in.system_size_y = 8;
    in.system_size_z = 8;
    in.lattice_type = SC;
    in.rand_seed = 0xDEADBEEF;
    in.overpotential_ramp_rate = 0.11;
    in.max_overpotential = 2.2;
    in.num_nn_levels = 1;
    in.num_bond_types = 1;
    in.num_elements = 1;
    in.atom_names_cnt = 0;
    in.nn_energy = NULL;
    in.dissolution = 0;

    initialize_env_from_input(&in, se);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(in.rand_seed, se->rand_seed, "rand_seed copied");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.11, se->overpotential_ramp_rate,
                                     "overpotential_ramp_rate copied");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(2.2, se->max_overpotential, "max_overpotential copied");
}

void test_initialize_simulation(void)
{
    struct UserInputs inputs = {0};
    inputs.system_size_x = 128;
    inputs.system_size_y = 128;
    inputs.system_size_z = 128;
    inputs.temperature = 293.;
    inputs.rand_seed = 12345;
    inputs.overpotential_ramp_rate = 0.03;
    inputs.max_overpotential = 1.2;
    inputs.num_nn_levels = 2;
    inputs.nn_energy = (double *)malloc(6 * sizeof(double));
    if (inputs.nn_energy == NULL) {
        perror("Couldn't allocate memory for nn_energy array");
        clean_and_error(errno);
    }
    inputs.nn_energy[0] = 0.15;
    inputs.nn_energy[1] = 0.15;
    inputs.nn_energy[2] = 0.15;
    inputs.nn_energy[3] = 0.10;
    inputs.nn_energy[4] = 0.10;
    inputs.nn_energy[5] = 0.10;

    ls->framenum = 0;
    inputs.lattice_type = FCC;
    inputs.geometry = GEOMETRY_CLUSTER;
    inputs.geometry_param = 32;
    inputs.num_elements = 2;
    inputs.num_bond_types = 3;
    inputs.atom_names_cnt = 2;
    inputs.atom_names = (char **)malloc(2 * sizeof(char *));
    if (inputs.atom_names == NULL) {
        perror("Couldn't allocate memory for atom_names array");
        clean_and_error(errno);
    }

    inputs.atom_names[0] = malloc(3);
    if (inputs.atom_names[0] == NULL) {
        perror("Couldn't allocate memory for atom_names array");
        clean_and_error(errno);
    }

    inputs.atom_names[1] = malloc(3);
    if (inputs.atom_names[1] == NULL) {
        perror("Couldn't allocate memory for atom_names array");
        clean_and_error(errno);
    }

    strcpy(inputs.atom_names[0], "Ag");
    strcpy(inputs.atom_names[1], "Au");

    inputs.dissolution = 1;

    inputs.is_soluble = (bool *)malloc(2 * sizeof(bool));
    if (inputs.is_soluble == NULL) {
        perror("Couldn't allocate memory for is_soluble array");
        clean_and_error(errno);
    }
    inputs.is_soluble[0] = 1;
    inputs.is_soluble[1] = 0;

    inputs.substrate_composition = (double *)malloc(2 * sizeof(double));
    if (inputs.substrate_composition == NULL) {
        perror("Couldn't allocate memory for substrate_composition array");
        clean_and_error(errno);
    }
    inputs.substrate_composition[0] = 0.75;
    inputs.substrate_composition[1] = 0.25;

    inputs.sim_end_type = SIM_END_BY_ITERATIONS;
    inputs.run_stime = 0;
    inputs.final_iteration = 2000;

    ls->sim_log = temp_log;

    initialize_simulation(&inputs, ss, se, ls);

    TEST_ASSERT_TRUE_MESSAGE(1, "initialize simulation");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293.0, ss->temperature, "ss.temperature set from inputs");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(12345u, se->rand_seed, "se.rand_seed copied from inputs");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.03, se->overpotential_ramp_rate,
                                     "se.overpotential_ramp_rate copied");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.2, se->max_overpotential, "se.max_overpotential copied");

    TEST_ASSERT_EQUAL_INT_MESSAGE(12, se->num_transition_vectors,
                                  "FCC should have 12 first-nearest transition vectors");
    TEST_ASSERT_EQUAL_INT_MESSAGE(18, se->num_energy_contributors,
                                  "With 2 nn levels FCC should have 18 energy contributors");

    TEST_ASSERT_TRUE_MESSAGE(se->max_atoms > 0, "se->max_atoms should be positive");

    TEST_ASSERT_TRUE_MESSAGE(ss->atom_cnt > 0, "initialize_initial_structure should add atoms");

    TEST_ASSERT_NOT_NULL_MESSAGE(ss->atom_arr, "ss->atom_arr should be allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->rate_arr, "ss->rate_arr should be allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->transition_arr, "ss->transition_arr should be allocated");

    TEST_ASSERT_NOT_NULL_MESSAGE(ss->transition_probability.rate_arr_index,
                                 "ss->transition_probability.rate_arr_index should be allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->transition_probability.lbound,
                                 "ss->transition_probability.lbound should be allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss->transition_probability.ubound,
                                 "ss->transition_probability.ubound should be allocated");

    TEST_ASSERT_NOT_NULL_MESSAGE(se->transition_vectors,
                                 "se->transition_vectors should be allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(se->opposite_tvectors,
                                 "se->opposite_tvectors should be allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(se->atoms_per_nn_level,
                                 "se->atoms_per_nn_level should be allocated");

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->is_soluble[0], "se->is_soluble[0] should be 1");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->is_soluble[1], "se->is_soluble[1] should be 0");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, se->substrate_composition[0],
                                     "substrate_composition[0]");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, se->substrate_composition[1],
                                     "substrate_composition[1]");
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_initialize_simulation);
    RUN_TEST(test_get_shifts_sets_expected_shifts);
    RUN_TEST(test_set_primitive_basis_fcc_sets_expected_vectors);
    RUN_TEST(test_initialize_neighbor_offsets_fcc_sets_expected_values);
    RUN_TEST(test_corners2limits_updates_limits_correctly);
    RUN_TEST(test_initialize_state_from_input_copies_values);
    RUN_TEST(test_initialize_zones_allocates_and_sets_offsets);
    RUN_TEST(test_initialize_simulation_box_and_lat_range);
    RUN_TEST(test_set_default_orientation_sc_returns_identity);
    RUN_TEST(test_allocate_simulation_arrays_allocates_expected_buffers);
    RUN_TEST(test_initialize_env_from_input_sets_seed_and_overpotentials);

    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
