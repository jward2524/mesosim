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

// mostly redundant with TestInput.c tests
void test_process_in_file_cluster(void)
{
    char filename[] = "test/cluster_nns.in";
    input_file = open_file(filename);

    process_in_file(input_file, ss, se, ls);
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
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->csv_framenum, "frame number");
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

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200., ls->csv_schedule.interval, "log interval");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200., ls->next_csv_checkpoint, "log checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->output_steps_csv, "output iter csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_state_csv, "output state csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_xyz, "output xyz");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_ITERATION, ls->csv_schedule.mode, "logging analysis type");
}

void test_process_in_file_mc(void)
{
    char filename[] = "test/mc.in";
    input_file = open_file(filename);

    process_in_file(input_file, ss, se, ls);
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
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_ITERATION, ls->csv_schedule.mode, "logging analysis type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->csv_framenum, "frame number");
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
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1., ls->csv_schedule.interval, "log interval");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1., ls->next_csv_checkpoint, "log checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->output_steps_csv, "output iter csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_state_csv, "output state csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->output_xyz, "output xyz");
}

void test_process_xyz_file(void)
{
    char filename[] = "test/sheet256.xyz";
    strncpy(se->atoms_filename, filename, strlen(filename));
    // initialize_from_file(ss, se, ls);

    atom_file = open_file(filename);

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

    process_xyz_file(atom_file, ss, se, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(499996, ss->atom_cnt, "Atom count");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293., ss->temperature, "Temperature");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.9, ss->overpotential, "Potential");
    TEST_ASSERT_EQUAL_INT_MESSAGE(190, ss->iter, "Iteration");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.3, ss->elapsed_stime, "Time");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(7.65, ss->total_internal_energy, "Energy");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, ls->xyz_framenum, "Frame");

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

    EXPECT_EXIT(1, {
        safe_log(temp_log, huge_format);
        TEST_FAIL_MESSAGE("Expected nnlevels argcount failure");
    });
}

// fflush failure (needs mocking, optional)
// Not trivial in standard C; usually requires dependency injection or linking fakes

void test_output_csv_header_success(void)
{
    ls->csv_field_count = 2;
    ls->csv_fields = malloc(2 * sizeof(char *));
    ls->csv_fields[0] = dup_str("iter");
    ls->csv_fields[1] = dup_str("energy");

    rewind(temp_log);
    output_csv_header(temp_log, ls);
    rewind(temp_log);

    char buffer[256];
    fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_EQUAL_STRING_MESSAGE("frame,iter,energy\n", buffer, "CSV header");
}

void test_log_state_csv_success(void)
{
    ss->iter = 42;
    ss->total_internal_energy = 3.14;
    
    ls->csv_framenum = 2;
    ls->csv_field_count = 2;
    ls->csv_field_funcs = malloc(2 * sizeof(CsvFieldFuncPtr));
    ls->csv_field_funcs[0] = csv_field_map[0].get_value; // get_iter
    ls->csv_field_funcs[1] = csv_field_map[2].get_value; // get_energy

    rewind(temp_log);
    log_state_csv(temp_log, ss, ls);
    rewind(temp_log);

    char buffer[256];
    fgets(buffer, sizeof(buffer), temp_log);
    // default precision of floats/doubles is 6
    TEST_ASSERT_EQUAL_STRING_MESSAGE("2,42,3.140000\n", buffer, "CSV state log");
}

void test_log_state_csv_one_field(void)
{
    ss->iter = 99;
    ls->csv_framenum = 5;
    ls->csv_field_count = 1;
    ls->csv_field_funcs = malloc(sizeof(CsvFieldFuncPtr));
    ls->csv_field_funcs[0] = csv_field_map[0].get_value; // get_iter

    rewind(temp_log);
    log_state_csv(temp_log, ss, ls);
    rewind(temp_log);

    char buffer[128];
    fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_EQUAL_STRING_MESSAGE("5,99\n", buffer, "CSV state log with one field");
}

void test_output_csv_header_one_field(void)
{
    ls->csv_field_count = 1;
    ls->csv_fields = malloc(sizeof(char *));
    ls->csv_fields[0] = dup_str("energy");

    rewind(temp_log);
    output_csv_header(temp_log, ls);
    rewind(temp_log);

    char buffer[128];
    fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_EQUAL_STRING_MESSAGE("frame,energy\n", buffer, "CSV header with one field");
}

void test_log_state_csv_mixed_fields(void)
{
    ss->iter = 7;
    ss->total_internal_energy = -1.23;
    ss->temperature = 273.15;
    ls->csv_framenum = 3;
    ls->csv_field_count = 3;
    ls->csv_field_funcs = malloc(3 * sizeof(CsvFieldFuncPtr));
    ls->csv_field_funcs[0] = csv_field_map[0].get_value; // get_iter
    ls->csv_field_funcs[1] = csv_field_map[2].get_value; // get_energy
    ls->csv_field_funcs[2] = csv_field_map[3].get_value; // get_temperature

    rewind(temp_log);
    log_state_csv(temp_log, ss, ls);
    rewind(temp_log);

    char buffer[256];
    fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_EQUAL_STRING_MESSAGE("3,7,-1.230000,273.150000\n", buffer, "CSV state log with mixed fields");
}

void test_output_csv_header_many_fields(void)
{
    size_t n = 10;
    ls->csv_field_count = n;
    ls->csv_fields = malloc((size_t)n * sizeof(char *));
    for (int i = 0; i < n; ++i) {
        char name[16];
        sprintf(name, "field%d", i);
        ls->csv_fields[i] = dup_str(name);
    }

    rewind(temp_log);
    output_csv_header(temp_log, ls);
    rewind(temp_log);

    char buffer[512];
    fgets(buffer, sizeof(buffer), temp_log);
    char expected[256] = "frame";
    for (int i = 0; i < n; ++i) {
        strcat(expected, ",");
        char name[16];
        sprintf(name, "field%d", i);
        strcat(expected, name);
    }
    strcat(expected, "\n");
    TEST_ASSERT_EQUAL_STRING_MESSAGE(expected, buffer, "CSV header with many fields");
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_process_in_file_mc);
    RUN_TEST(test_process_in_file_cluster);
    RUN_TEST(test_process_xyz_file);
    RUN_TEST(test_safe_log_writes);
    RUN_TEST(test_safe_log_long_line);
    RUN_TEST(test_safe_log_buffer_overflow);
    
    RUN_TEST(test_output_csv_header_success);
    RUN_TEST(test_log_state_csv_success);
    RUN_TEST(test_output_csv_header_one_field);
    RUN_TEST(test_output_csv_header_many_fields);
    RUN_TEST(test_log_state_csv_one_field);
    RUN_TEST(test_log_state_csv_mixed_fields);

    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
