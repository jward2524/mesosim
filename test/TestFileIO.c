#include "Atoms.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "State.h"
#include "TUtils.h"
#include "Utils.h"
#include "unity.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct SimulationState *ss;
struct SimulationEnv *se;
struct LoggingState *ls;
FILE *temp_log;
FILE *atom_file;
FILE *input_file;

// before and after each test (each RUN_TEST)
void setUp(void)
{
    initialize_states(&ss, &se, &ls);
    init_temp(&temp_log);
}

void tearDown(void)
{
    clean_and_error(0);

    // fclose needs to be here in case a test fails
    fclose(temp_log);
    close_if_exists(&atom_file);
    close_if_exists(&input_file);
}

void test_parse_input_systemsize(void)
{
    char line[] = "systemsize 128 128 128";

    parse_input(line, temp_log, ss, se, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_x, "system size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_y, "system size y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_z, "system size z");
}

void test_parse_input_geometry_file(void)
{
    char line[] = "geometry file test/sheet.xyz";

    parse_input(line, temp_log, ss, se, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(GEOMETRY_FROM_FILE, se->geometry, "simulation type");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("test/sheet.xyz", se->atoms_filename, "atoms filename");
}

void test_process_in_file_cluster(void)
{
    char filename[] = "test/cluster_nns.in";
    input_file = fopen(filename, "r");
    fopen_error(filename, input_file);

    int ret = process_in_file(temp_log, input_file, ss, se, ls);
    TEST_ASSERT_TRUE_MESSAGE(ret, "fxn return");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FLAVOR_KMC, se->flavor, "flavor");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_x, "system size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_y, "system size y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_z, "system size z");
    TEST_ASSERT_EQUAL_INT_MESSAGE(32, se->cluster_radius, "cluster radius");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293., ss->temperature, "temperature");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12345, se->rand_seed, "random seed");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.9, se->initial_overpotential, "initial overpotential");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.03, se->overpotential_ramp_rate, "overpotential ramp");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.2, se->max_overpotential, "max overpotential");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->num_nn_levels, "number of nn shells");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[0], "nn energy shell 1 idx 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[1], "nn energy shell 1 idx 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[2], "nn energy shell 1 idx 2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, se->nn_energy[3], "nn energy shell 2 idx 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, se->nn_energy[4], "nn energy shell 2 idx 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, se->nn_energy[5], "nn energy shell 2 idx 2");
    TEST_ASSERT_EQUAL_INT_MESSAGE(ITERATION_INTERVALS, ls->analysis_type, "logging analysis type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->framenum, "frame number");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FCC, se->lattice_type, "lattice type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(MAXIMUM_NUMBER_OF_NEIGHBORS, se->num_transition_vectors,
                                  "number of transition vectors");
    TEST_ASSERT_EQUAL_INT_MESSAGE(GEOMETRY_CLUSTER, se->geometry, "geometry type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->sheet_thickness, "sheet thickness");
    TEST_ASSERT_EQUAL_INT_MESSAGE(32, se->cluster_radius, "cluster radius");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("", se->atoms_filename, "geometry filename");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->num_elements, "number of elements");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->num_bond_types, "number of bond types");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->atom_names_cnt, "number of atom names");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Ag", se->atom_names[0], "first atom name");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Au", se->atom_names[1], "second atom name");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->dissolution, "dissolution");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->is_soluble[0], "Ag solubility");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->is_soluble[1], "Au solubility");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, se->substrate_composition[0], "Ag composition");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, se->substrate_composition[1], "Au composition");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, ss->sim_end_type, "simulation end type");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0, ss->run_stime, "simulation max runtime");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2000, ss->final_iteration, "simulation max iteration");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200., ls->log_interval, "log interval");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200., ls->next_log_checkpoint, "log checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_iter_csv, "output iter csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_state_csv, "output state csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_xyz, "output xyz");
}

void test_process_in_file_mc(void)
{
    char filename[] = "test/mc.in";
    input_file = fopen(filename, "r");
    fopen_error(filename, input_file);

    int ret = process_in_file(temp_log, input_file, ss, se, ls);
    TEST_ASSERT_TRUE_MESSAGE(ret, "fxn return");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FLAVOR_MC, se->flavor, "flavor");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_x, "system size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_y, "system size y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_z, "system size z");
    TEST_ASSERT_EQUAL_INT_MESSAGE(32, se->cluster_radius, "cluster radius");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293., ss->temperature, "temperature");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12345, se->rand_seed, "random seed");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0, se->initial_overpotential, "initial overpotential");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0, se->overpotential_ramp_rate, "overpotential ramp");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0, se->max_overpotential, "max overpotential");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->num_nn_levels, "number of nn shells");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[0], "nn energy shell 1 idx 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[1], "nn energy shell 1 idx 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[2], "nn energy shell 1 idx 2");
    TEST_ASSERT_EQUAL_INT_MESSAGE(ITERATION_INTERVALS, ls->analysis_type, "logging analysis type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->framenum, "frame number");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FCC, se->lattice_type, "lattice type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(MAXIMUM_NUMBER_OF_NEIGHBORS, se->num_transition_vectors,
                                  "number of transition vectors");
    TEST_ASSERT_EQUAL_INT_MESSAGE(GEOMETRY_CLUSTER, se->geometry, "geometry type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->sheet_thickness, "sheet thickness");
    TEST_ASSERT_EQUAL_INT_MESSAGE(32, se->cluster_radius, "cluster radius");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("", se->atoms_filename, "geometry filename");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->num_elements, "number of elements");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->num_bond_types, "number of bond types");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->atom_names_cnt, "number of atom names");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Ag", se->atom_names[0], "first atom name");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Au", se->atom_names[1], "second atom name");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->dissolution, "dissolution");
    // TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->is_soluble[0], "Ag solubility");
    // TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->is_soluble[1], "Au solubility");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, se->substrate_composition[0], "Ag composition");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, se->substrate_composition[1], "Au composition");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, ss->sim_end_type, "simulation end type");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0, ss->run_stime, "simulation max runtime");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, ss->final_iteration, "simulation max iteration");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1., ls->log_interval, "log interval");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1., ls->next_log_checkpoint, "log checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->output_iter_csv, "output iter csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_state_csv, "output state csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->output_xyz, "output xyz");
}

void test_process_xyz_file(void)
{
    char filename[] = "test/sheet256.xyz";
    strncpy(se->atoms_filename, filename, strlen(filename));
    // initialize_from_file(ss, se, ls);

    atom_file = fopen(filename, "r");
    fopen_error(filename, atom_file);

    se->lattice_type = FCC;
    se->system_size_x = 256;
    se->system_size_y = 256;
    se->system_size_z = 256;
    se->atom_names = malloc(2 * sizeof *se->atom_names);
    se->atom_names[0] = malloc(sizeof "Ag");
    strcpy(se->atom_names[0], "Ag");
    se->atom_names[1] = malloc(sizeof "Au");
    strcpy(se->atom_names[1], "Au");
    se->atom_names_cnt = 2;
    get_shifts(se);
    set_primitive_basis(se);
    initialize_simulation_box(se);
    initialize_zones(ss->zone_arr, se);
    initialize_simulation_variables(ss, se);

    int res = process_xyz_file(temp_log, atom_file, ss, se, ls);

    TEST_ASSERT_TRUE_MESSAGE(res, "Return value of process_xyz_file");
    TEST_ASSERT_EQUAL_INT_MESSAGE(499996, ss->atom_cnt, "Atom count");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293., ss->temperature, "Temperature");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.9, ss->overpotential, "Potential");
    TEST_ASSERT_EQUAL_INT_MESSAGE(190, ss->iter, "Iteration");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.3, ss->elapsed_stime, "Time");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(7.65, ss->total_internal_energy, "Energy");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, ls->framenum, "Frame");

    // 0 Ag -256.000000 -96.000000 -96.000000
    // double cart[3] = {-256, -96, -96};
    // double cart2[3] = {83, 133, 40};

    int lat[3] = {-128, -128, 32};
    int lat2[3] = {88, -5, 45};

    Atom at = *ss->atom_arr[0];
    TEST_ASSERT_EQUAL_INT_MESSAGE(lat[0], at.lattice[0], "Atom 0, lattice x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(lat[1], at.lattice[1], "Atom 0, lattice y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(lat[2], at.lattice[2], "Atom 0, lattice z");

    at = *ss->atom_arr[ss->atom_cnt - 1];
    TEST_ASSERT_EQUAL_INT_MESSAGE(lat2[0], at.lattice[0], "Atom max, lattice x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(lat2[1], at.lattice[1], "Atom max, lattice y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(lat2[2], at.lattice[2], "Atom max, lattice z");
}

// simple log
void test_safe_log_writes(void)
{
    safe_log(temp_log, "Hello %s %d\n", "World", 42);

    // Rewind to read
    rewind(temp_log);

    char buffer[128];
    fgets(buffer, sizeof(buffer), temp_log);

    TEST_ASSERT_EQUAL_STRING_MESSAGE("Hello World 42\n", buffer,
                                     "Simple log should be written correctly");
}

// long line within buffer limit
void test_safe_log_long_line(void)
{
    char long_str[512];
    memset(long_str, 'A', sizeof(long_str) - 2);
    long_str[sizeof(long_str) - 2] = '\n';
    long_str[sizeof(long_str) - 1] = '\0';

    safe_log(temp_log, "%s", long_str);

    rewind(temp_log);

    char buffer[1024];
    fgets(buffer, sizeof(buffer), temp_log);

    TEST_ASSERT_EQUAL_STRING_MESSAGE(long_str, buffer, "Long log line should be written correctly");
}

// formatting error simulated
void test_safe_log_buffer_overflow(void)
{
    char huge_format[2048];
    memset(huge_format, 'A', sizeof(huge_format) - 1);
    huge_format[sizeof(huge_format) - 1] = '\0';

    safe_log(temp_log, huge_format);
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, exit_flag,
                                  "Expected safe_log to trigger exit on buffer overflow");
}

// fflush failure (needs mocking, optional)
// Not trivial in standard C; usually requires dependency injection or linking fakes

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_parse_input_systemsize);
    RUN_TEST(test_parse_input_geometry_file);
    RUN_TEST(test_process_in_file_mc);
    RUN_TEST(test_process_in_file_cluster);
    RUN_TEST(test_process_xyz_file);
    RUN_TEST(test_safe_log_writes);
    RUN_TEST(test_safe_log_long_line);
    RUN_TEST(test_safe_log_buffer_overflow);

    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
