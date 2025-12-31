#include "Initialization.h"
#include "unity.h"
#include <stdlib.h>

struct SimulationState *ss;
struct SimulationEnv *se;
struct LoggingState *ls;

void setUp(void)
{
    ss = calloc(1, sizeof(struct SimulationState));
    se = calloc(1, sizeof(struct SimulationEnv));
    ls = calloc(1, sizeof(struct LoggingState));
}

void tearDown(void)
{
    free(se->nn_energy);
    free(se->atom_names);
    free(se->is_soluble);
    free(se->substrate_composition);
    free(ss);
    free(se);
    free(ls);
}

void test_initialize_simulation(void)
{
    se->system_size_x = 128;
    se->system_size_y = 128;
    se->system_size_z = 128;
    se->cluster_radius = 32;
    ss->temperature = 293.;
    se->rand_seed = 12345;
    se->initial_overpotential = 0.9;
    se->overpotential_ramp_rate = 0.03;
    se->max_overpotential = 1.2;
    se->num_nn_levels = 2;
    se->nn_energy = (double *)malloc(6 * sizeof(double));
    se->nn_energy[0] = 0.15;
    se->nn_energy[1] = 0.15;
    se->nn_energy[2] = 0.15;
    se->nn_energy[3] = 0.10;
    se->nn_energy[4] = 0.10;
    se->nn_energy[5] = 0.10;
    ls->analysis_type = ITERATION_INTERVALS;
    ls->framenum = 0;
    se->lattice_type = FCC;
    se->num_transition_vectors = MAXIMUM_NUMBER_OF_NEIGHBORS;
    se->geometry = GEOMETRY_CLUSTER;
    se->sheet_thickness = 0;
    se->cluster_radius = 32;
    // se->atoms_filename = "";
    se->num_elements = 2;
    se->num_bond_types = 3;
    se->atom_names_cnt = 2;
    se->atom_names = (char **)malloc(2 * sizeof(char*));
    se->atom_names[0] = "Ag";
    se->atom_names[1] = "Au";
    se->dissolution = 1;
    se->is_soluble = (bool *)malloc(2 * sizeof(bool));
    se->is_soluble[0] = 1;
    se->is_soluble[1] = 0;
    se->substrate_composition = (double *)malloc(2 * sizeof(double));
    se->substrate_composition[0] = 0.75;
    se->substrate_composition[1] = 0.25;
    ss->sim_end_type = SIM_END_BY_ITERATIONS;
    ss->run_stime = 0;
    ss->final_iteration = 2000;
    
    initialize_simulation(ss, se, ls);
    TEST_ASSERT_TRUE_MESSAGE(1, "initialize simulation");
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_initialize_simulation);
    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
