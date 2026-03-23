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
    if (se->nn_energy == NULL) {
        perror("Couldn't allocate memory for nn_energy array");
        clean_and_error(errno);
    }
    se->nn_energy[0] = 0.15;
    se->nn_energy[1] = 0.15;
    se->nn_energy[2] = 0.15;
    se->nn_energy[3] = 0.10;
    se->nn_energy[4] = 0.10;
    se->nn_energy[5] = 0.10;

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
    se->atom_names = (char **)malloc(2 * sizeof(char *));
    if (se->atom_names == NULL) {
        perror("Couldn't allocate memory for atom_names array");
        clean_and_error(errno);
    }

    se->atom_names[0] = malloc(3);
    if (se->atom_names[0] == NULL) {
        perror("Couldn't allocate memory for atom_names array");
        clean_and_error(errno);
    }

    se->atom_names[1] = malloc(3);
    if (se->atom_names[1] == NULL) {
        perror("Couldn't allocate memory for atom_names array");
        clean_and_error(errno);
    }

    strcpy(se->atom_names[0], "Ag");
    strcpy(se->atom_names[1], "Au");

    se->dissolution = 1;

    se->is_soluble = (bool *)malloc(2 * sizeof(bool));
    if (se->is_soluble == NULL) {
        perror("Couldn't allocate memory for is_soluble array");
        clean_and_error(errno);
    }
    se->is_soluble[0] = 1;
    se->is_soluble[1] = 0;

    se->substrate_composition = (double *)malloc(2 * sizeof(double));
    if (se->substrate_composition == NULL) {
        perror("Couldn't allocate memory for substrate_composition array");
        clean_and_error(errno);
    }
    se->substrate_composition[0] = 0.75;
    se->substrate_composition[1] = 0.25;

    ss->sim_end_type = SIM_END_BY_ITERATIONS;
    ss->run_stime = 0;
    ss->final_iteration = 2000;

    ls->sim_log = temp_log;

    initialize_simulation(ss, se, ls);
    TEST_ASSERT_TRUE_MESSAGE(1, "initialize simulation");
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_initialize_simulation);
    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
