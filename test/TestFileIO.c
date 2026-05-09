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
    // fclose needs to be here in case a test fails
    // set to null to prevent double free from fclose + clean_and_error
    if (gp_log_state && (*gp_log_state)->sim_log) {
        (*gp_log_state)->sim_log = NULL;
    }
    fclose(temp_log);
    clean_and_error(0);
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

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_state_csv, "output state csv");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200., ls->csv_schedule.interval, "csv log interval");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200., ls->next_csv_checkpoint, "csv log checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->output_steps_csv, "output steps csv");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_ITERATION, ls->csv_schedule.mode,
                                  "logging analysis type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->output_xyz, "output xyz");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->xyz_stripped, "stripped xyz");
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
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_ITERATION, ls->csv_schedule.mode,
                                  "logging analysis type");
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
    strncpy(se->atoms_filename, filename, strlen(se->atoms_filename));
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
    initialize_zones(&ss->zone_arr, se);
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
    char *ret = fgets(buffer, sizeof(buffer), temp_log);

    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
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
    char *ret = fgets(buffer, sizeof(buffer), temp_log);

    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
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

// get_precision tests - core functionality
void test_get_precision_positive_numbers(void)
{
    // total=100, increment=1, incr_precision=2 -> log10(100)=2, log10(1)=0, diff=2, result=4
    int result = get_precision(100.0, 1.0, 2);
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, result, "positive numbers");
}

void test_get_precision_large_difference(void)
{
    // total=1000000, increment=1, incr_precision=1 -> log10(1e6)=6, log10(1)=0, diff=6, result=7
    int result = get_precision(1000000.0, 1.0, 1);
    TEST_ASSERT_EQUAL_INT_MESSAGE(7, result, "large magnitude difference");
}

void test_get_precision_negative_total(void)
{
    // total=-100, increment=1, incr_precision=2 -> uses abs value, result=4
    int result = get_precision(-100.0, 1.0, 2);
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, result, "negative total");
}

void test_get_precision_negative_increment(void)
{
    // total=100, increment=-1, incr_precision=2 -> uses abs value, result=4
    int result = get_precision(100.0, -1.0, 2);
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, result, "negative increment");
}

void test_get_precision_both_negative(void)
{
    // total=-50, increment=-0.5, incr_precision=1 -> log10(50)=1, log10(0.5)=-0, diff=1, result=2
    int result = get_precision(-50.0, -0.5, 1);
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, result, "both negative");
}

void test_get_precision_total_zero(void)
{
    // When total is zero, should return incr_precision without computing logs
    int result = get_precision(0.0, 5.0, 3);
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, result, "total is zero");
}

void test_get_precision_increment_zero(void)
{
    // When increment is zero, should return incr_precision without computing logs
    int result = get_precision(100.0, 0.0, 4);
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, result, "increment is zero");
}

void test_get_precision_zero_incr_precision(void)
{
    // total=100, increment=1, incr_precision=0 -> log diff=2, result=2
    int result = get_precision(100.0, 1.0, 0);
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, result, "incr_precision is zero");
}

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
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
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
    ls->csv_fields = malloc(2 * sizeof(char *));
    ls->csv_fields[0] = dup_str("iter");
    ls->csv_fields[1] = dup_str("energy");

    rewind(temp_log);
    log_state_csv(temp_log, ss, ls);
    rewind(temp_log);

    char buffer[256];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
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
    ls->csv_fields = malloc(1 * sizeof(char *));
    ls->csv_fields[0] = dup_str("iter");

    rewind(temp_log);
    log_state_csv(temp_log, ss, ls);
    rewind(temp_log);

    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
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
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("frame,energy\n", buffer, "CSV header with one field");
}

void test_log_state_csv_mixed_fields(void)
{
    ss->iter = 7;
    ss->total_internal_energy = -1.23;
    ss->overpotential = 0.9;
    ls->overpot_precision = 4;
    ls->csv_framenum = 3;
    ls->csv_field_count = 3;
    ls->csv_field_funcs = malloc(3 * sizeof(CsvFieldFuncPtr));
    ls->csv_field_funcs[0] = csv_field_map[0].get_value; // get_iter
    ls->csv_field_funcs[1] = csv_field_map[2].get_value; // get_energy
    ls->csv_field_funcs[2] = csv_field_map[4].get_value; // get_overpotential
    ls->csv_fields = malloc(3 * sizeof(char *));
    ls->csv_fields[0] = dup_str("iter");
    ls->csv_fields[1] = dup_str("energy");
    ls->csv_fields[2] = dup_str("overpotential");

    rewind(temp_log);
    log_state_csv(temp_log, ss, ls);
    rewind(temp_log);

    char buffer[256];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("3,7,-1.230000,9.0000e-01\n", buffer,
                                     "CSV state log with mixed fields");
}

void test_output_csv_header_many_fields(void)
{
    int n = 10;
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
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");

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

void test_output_mc_steps_header_coord(void)
{
    bool output_coord = true;
    rewind(temp_log);
    output_mc_steps_header(temp_log, output_coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE(
        "iter,energy,deltaE,performed,u1,v1,w1,u2,v2,w2,coordination\n", buffer,
        "MC iter header coord");
}

void test_output_mc_steps_header_coordless(void)
{
    bool output_coord = false;
    rewind(temp_log);
    output_mc_steps_header(temp_log, output_coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("iter,energy,deltaE,performed,u1,v1,w1,u2,v2,w2\n", buffer,
                                     "MC iter header coordless");
}

void test_log_mc_steps_basic_coordless(void)
{
    unsigned long int iter = 123;
    double sys_energy = 4.56;
    double deltaE = -0.12;
    int performed = 1;
    int uvw1[3] = {1, 2, 3};
    int uvw2[3] = {4, 5, 6};
    int coord = -1;

    rewind(temp_log);
    log_mc_steps(temp_log, iter, sys_energy, deltaE, performed, uvw1, uvw2, coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("123,4.560000,-0.120000,1,1,2,3,4,5,6\n", buffer,
                                     "MC iter log basic coordless");
}

void test_log_mc_steps_basic_coord(void)
{
    unsigned long int iter = 123;
    double sys_energy = 4.56;
    double deltaE = -0.12;
    int performed = 1;
    int uvw1[3] = {1, 2, 3};
    int uvw2[3] = {4, 5, 6};
    int coord = 9;

    rewind(temp_log);
    log_mc_steps(temp_log, iter, sys_energy, deltaE, performed, uvw1, uvw2, coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("123,4.560000,-0.120000,1,1,2,3,4,5,6,9\n", buffer,
                                     "MC iter log basic coord");
}

void test_output_kmc_steps_header_coord(void)
{
    bool output_coord = true;
    rewind(temp_log);
    output_kmc_steps_header(temp_log, output_coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("iter,sim_time,energy,u1,v1,w1,u2,v2,w2,coordination\n",
                                     buffer, "KMC iter header coord");
}

void test_output_kmc_steps_header_coordless(void)
{
    bool output_coord = false;
    rewind(temp_log);
    output_kmc_steps_header(temp_log, output_coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("iter,sim_time,energy,u1,v1,w1,u2,v2,w2\n", buffer,
                                     "KMC iter header coordless");
}

void test_log_kmc_steps_basic_coord(void)
{
    unsigned long int iter = 42;
    double sim_time = 0.5;
    double sys_energy = 7.89;
    int uvw1[3] = {7, 8, 9};
    int uvw2[3] = {10, 11, 12};
    int is_evap = 0;
    int coord = 2;

    rewind(temp_log);
    log_kmc_steps(temp_log, iter, sim_time, 6, sys_energy, uvw1, uvw2, is_evap, coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("42,5.000000e-01,7.890000,7,8,9,10,11,12,2\n", buffer,
                                     "KMC iter log basic coord");
}

void test_log_kmc_steps_evap(void)
{
    unsigned long int iter = 99;
    double sim_time = 1.23;
    double sys_energy = -2.34;
    int uvw1[3] = {3, 2, 1};
    int uvw2[3] = {0, 0, 0};
    int is_evap = 1;
    int coord = -1;

    rewind(temp_log);
    log_kmc_steps(temp_log, iter, sim_time, 6, sys_energy, uvw1, uvw2, is_evap, coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("99,1.230000e+00,-2.340000,3,2,1,,,\n", buffer,
                                     "KMC iter log evaporation");
}

void test_log_kmc_steps_evap_coord(void)
{
    unsigned long int iter = 99;
    double sim_time = 1.23;
    double sys_energy = -2.34;
    int uvw1[3] = {3, 2, 1};
    int uvw2[3] = {0, 0, 0};
    int is_evap = 1;
    int coord = 8;

    rewind(temp_log);
    log_kmc_steps(temp_log, iter, sim_time, 6, sys_energy, uvw1, uvw2, is_evap, coord);
    rewind(temp_log);
    char buffer[128];
    char *ret = fgets(buffer, sizeof(buffer), temp_log);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("99,1.230000e+00,-2.340000,3,2,1,,,,8\n", buffer,
                                     "KMC iter log evaporation");
}

void test_write_xyz_suffix_iteration(void)
{
    char suffix[256];
    write_xyz_suffix(suffix, OUTPUT_SCHEDULE_INTERVAL_ITERATION, 1234);
    TEST_ASSERT_EQUAL_STRING_MESSAGE("i1234", suffix, "Suffix for iteration mode");
    write_xyz_suffix(suffix, OUTPUT_SCHEDULE_LIST_ITERATION, 5678);
    TEST_ASSERT_EQUAL_STRING_MESSAGE("i5678", suffix, "Suffix for list iteration mode");
}

void test_write_xyz_suffix_time(void)
{
    char suffix[256];
    write_xyz_suffix(suffix, OUTPUT_SCHEDULE_INTERVAL_TIME, 1.2345);
    TEST_ASSERT_EQUAL_STRING_MESSAGE("t1.2345", suffix, "Suffix for time mode");
    write_xyz_suffix(suffix, OUTPUT_SCHEDULE_LIST_TIME, 0.0001);
    TEST_ASSERT_EQUAL_STRING_MESSAGE("t0.0001", suffix, "Suffix for list time mode");
}

void test_fstring_to_buffer_returns_null_on_format_error(void)
{
    // Intentionally pass an invalid format string
    // errors depend on implementation
    const char *result = fstring_to_buffer("%q", 42);
    TEST_ASSERT_NULL_MESSAGE(result, "Should return NULL on format error");
}

void test_input_logging_basic(void)
{
    se->system_size_x = 10;
    se->system_size_y = 20;
    se->system_size_z = 30;
    se->lattice_type = FCC;
    se->num_elements = 1;
    se->atom_names = malloc(sizeof(char *));
    se->atom_names[0] = dup_str("Ag");
    se->substrate_composition = malloc(sizeof(double));
    se->substrate_composition[0] = 1.0;
    se->is_soluble = malloc(sizeof(int));
    se->is_soluble[0] = 1;
    se->num_nn_levels = 1;
    se->nn_energy = malloc(sizeof(double));
    se->nn_energy[0] = 0.1;
    se->geometry = GEOMETRY_CLUSTER;
    se->cluster_radius = 5;
    se->rand_seed = 123;
    ss->atom_cnt = 42;
    ss->temperature = 300.0;
    ss->overpotential = 0.5;
    ls->sim_log = temp_log;
    ls->output_state_csv = 0;
    ls->output_xyz = 0;

    rewind(temp_log);
    input_logging(ss, se, ls);
    rewind(temp_log);

    ss->atom_cnt = 0;        // prevent free of uninitialized atoms in tearDown
    free(se->atom_names[0]); // free dup_str malloc'd string

    char buffer[1024];
    size_t bufstr_len = sizeof(buffer) - 1;
    size_t ret = fread(buffer, 1, bufstr_len, temp_log);
    TEST_ASSERT_TRUE_MESSAGE(ret > 0, "fread should have read data");
    buffer[bufstr_len] = '\0';

    // check that info is a substring of log output
    TEST_ASSERT_NOT_NULL(strstr(buffer, "System size is 10 x 20 x 30"));
    TEST_ASSERT_NOT_NULL(strstr(buffer, "Crystal structure is FCC"));
    TEST_ASSERT_NOT_NULL(strstr(buffer, "Random seed is 123"));
    TEST_ASSERT_NOT_NULL(strstr(buffer, "Initialized spherical cluster with radius 5"));
    TEST_ASSERT_NOT_NULL(strstr(buffer, "Atoms created, 42 total"));
}

// TODO: write_xyz and write_logs tests need simulation variables to be initialized
void test_write_xyz_file_creates_file_and_content(void)
{
    char prefix[] = "test_xyz_output";
    char suffix[] = "mysuffix";
    int frame_num = 1;
    int stripped = 0;
    ss->atom_cnt = 1;
    ss->iter = 2;
    ss->elapsed_stime = 3.0;
    ss->temperature = 4.0;
    ss->overpotential = 5.0;
    ss->total_internal_energy = 6.0;
    ss->atom_arr = malloc(sizeof(Atom *));
    ss->atom_arr[0] = malloc(sizeof(Atom));
    ss->atom_arr[0]->type = 0;
    ss->atom_arr[0]->cartesian[0] = 1.1;
    ss->atom_arr[0]->cartesian[1] = 2.2;
    ss->atom_arr[0]->cartesian[2] = 3.3;
    se->atom_names = malloc(sizeof(char *));
    se->atom_names[0] = dup_str("Ag");
    se->num_transition_vectors = 1;
    se->simbox_vectors_cart[0][0] = 1.0;
    se->simbox_vectors_cart[0][1] = 0.0;
    se->simbox_vectors_cart[0][2] = 0.0;
    se->simbox_vectors_cart[1][0] = 0.0;
    se->simbox_vectors_cart[1][1] = 1.0;
    se->simbox_vectors_cart[1][2] = 0.0;
    se->simbox_vectors_cart[2][0] = 0.0;
    se->simbox_vectors_cart[2][1] = 0.0;
    se->simbox_vectors_cart[2][2] = 1.0;
    se->simbox_origin_cart[0] = 0.0;
    se->simbox_origin_cart[1] = 0.0;
    se->simbox_origin_cart[2] = 0.0;

    bool result = write_xyz_file(prefix, frame_num, suffix, stripped, ss, se);

    TEST_ASSERT_TRUE(result);
    char filename[520];
    sprintf(filename, "%s_%d_%s.xyz", prefix, frame_num, suffix);
    FILE *f = fopen(filename, "r");
    TEST_ASSERT_NOT_NULL(f);
    char line[256];
    char *ret = fgets(line, sizeof(line), f);
    TEST_ASSERT_NOT_NULL_MESSAGE(ret, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING("1\n", line);
    // Clean up
    fclose(f);
    remove(filename);
}

void test_write_logs_increments_framenums(void)
{
    ls->sim_log = temp_log;
    ls->state_csv = temp_log;
    strncpy(ls->xyz_prefix, "test_xyz", sizeof(ls->xyz_prefix));
    strncpy(ls->xyz_suffix, "mysuffix", sizeof(ls->xyz_suffix));
    ls->xyz_stripped = 0;
    ls->csv_framenum = 0;
    ls->xyz_framenum = 0;
    ls->framenum = 0;
    ss->atom_cnt = 1;
    ss->iter = 1;
    ss->elapsed_stime = 1.0;
    ss->temperature = 1.0;
    ss->overpotential = 1.0;
    ss->total_internal_energy = 1.0;
    ss->atom_arr = malloc(sizeof(Atom *));
    ss->atom_arr[0] = malloc(sizeof(Atom));
    ss->atom_arr[0]->type = 0;
    ss->atom_arr[0]->cartesian[0] = 0.0;
    ss->atom_arr[0]->cartesian[1] = 0.0;
    ss->atom_arr[0]->cartesian[2] = 0.0;
    se->atom_names = malloc(sizeof(char *));
    se->atom_names[0] = dup_str("Ag");
    se->num_transition_vectors = 1;
    se->simbox_vectors_cart[0][0] = 1.0;
    se->simbox_vectors_cart[0][1] = 0.0;
    se->simbox_vectors_cart[0][2] = 0.0;
    se->simbox_vectors_cart[1][0] = 0.0;
    se->simbox_vectors_cart[1][1] = 1.0;
    se->simbox_vectors_cart[1][2] = 0.0;
    se->simbox_vectors_cart[2][0] = 0.0;
    se->simbox_vectors_cart[2][1] = 0.0;
    se->simbox_vectors_cart[2][2] = 1.0;
    se->simbox_origin_cart[0] = 0.0;
    se->simbox_origin_cart[1] = 0.0;
    se->simbox_origin_cart[2] = 0.0;

    write_logs(1, 1, ss, se, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->framenum, "frame number");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->csv_framenum, "csv frame number");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->xyz_framenum, "xyz frame number");
    // Clean up xyz file
    char filename[520];
    sprintf(filename, "%s_%d_%s.xyz", ls->xyz_prefix, 0, ls->xyz_suffix);
    remove(filename);
}

void test_output_if_passed_checkpoint_triggers_write_logs(void)
{
    ls->sim_log = temp_log;
    ls->state_csv = temp_log;
    strncpy(ls->xyz_prefix, "test_xyz", sizeof(ls->xyz_prefix));
    strncpy(ls->xyz_suffix, "mysuffix", sizeof(ls->xyz_suffix));
    ls->xyz_stripped = 0;
    ls->csv_framenum = 0;
    ls->xyz_framenum = 0;
    ls->framenum = 0;
    ss->atom_cnt = 1;
    ss->iter = 10;
    ss->elapsed_stime = 1.0;
    ss->temperature = 1.0;
    ss->overpotential = 1.0;
    ss->total_internal_energy = 1.0;
    ss->atom_arr = malloc(sizeof(Atom *));
    ss->atom_arr[0] = malloc(sizeof(Atom));
    ss->atom_arr[0]->type = 0;
    ss->atom_arr[0]->cartesian[0] = 0.0;
    ss->atom_arr[0]->cartesian[1] = 0.0;
    ss->atom_arr[0]->cartesian[2] = 0.0;
    se->atom_names = malloc(sizeof(char *));
    se->atom_names[0] = dup_str("Ag");
    se->num_transition_vectors = 1;
    se->simbox_vectors_cart[0][0] = 1.0;
    se->simbox_vectors_cart[0][1] = 0.0;
    se->simbox_vectors_cart[0][2] = 0.0;
    se->simbox_vectors_cart[1][0] = 0.0;
    se->simbox_vectors_cart[1][1] = 1.0;
    se->simbox_vectors_cart[1][2] = 0.0;
    se->simbox_vectors_cart[2][0] = 0.0;
    se->simbox_vectors_cart[2][1] = 0.0;
    se->simbox_vectors_cart[2][2] = 1.0;
    se->simbox_origin_cart[0] = 0.0;
    se->simbox_origin_cart[1] = 0.0;
    se->simbox_origin_cart[2] = 0.0;
    ls->output_state_csv = 1;
    ls->output_xyz = 1;
    ls->csv_schedule.mode = OUTPUT_SCHEDULE_LIST_ITERATION;
    double csv_list[1] = {10};
    ls->csv_schedule.list = csv_list;
    ls->csv_schedule.list_len = 1;
    ls->csv_schedule.list_idx = 0;
    ls->next_csv_checkpoint = 10;
    ls->xyz_schedule.mode = OUTPUT_SCHEDULE_LIST_ITERATION;
    double xyz_list[1] = {10};
    ls->xyz_schedule.list = xyz_list;
    ls->xyz_schedule.list_len = 1;
    ls->xyz_schedule.list_idx = 0;
    ls->next_xyz_checkpoint = 10;

    output_if_passed_checkpoint(ss, se, ls);

    TEST_ASSERT_EQUAL_INT(1, ls->framenum);
    TEST_ASSERT_EQUAL_INT(1, ls->csv_framenum);
    TEST_ASSERT_EQUAL_INT(1, ls->xyz_framenum);
    // Clean up xyz file
    char filename[520];
    sprintf(filename, "%s_%d_%s.xyz", ls->xyz_prefix, 0, ls->xyz_suffix);
    remove(filename);
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

    RUN_TEST(test_get_precision_positive_numbers);
    RUN_TEST(test_get_precision_large_difference);
    RUN_TEST(test_get_precision_negative_total);
    RUN_TEST(test_get_precision_negative_increment);
    RUN_TEST(test_get_precision_both_negative);
    RUN_TEST(test_get_precision_total_zero);
    RUN_TEST(test_get_precision_increment_zero);
    RUN_TEST(test_get_precision_zero_incr_precision);

    RUN_TEST(test_output_csv_header_success);
    RUN_TEST(test_log_state_csv_success);
    RUN_TEST(test_output_csv_header_one_field);
    RUN_TEST(test_output_csv_header_many_fields);
    RUN_TEST(test_log_state_csv_one_field);
    RUN_TEST(test_log_state_csv_mixed_fields);

    RUN_TEST(test_output_mc_steps_header_coord);
    RUN_TEST(test_output_mc_steps_header_coordless);
    RUN_TEST(test_log_mc_steps_basic_coordless);
    RUN_TEST(test_log_mc_steps_basic_coord);
    RUN_TEST(test_output_kmc_steps_header_coord);
    RUN_TEST(test_output_kmc_steps_header_coordless);
    RUN_TEST(test_log_kmc_steps_basic_coord);
    RUN_TEST(test_log_kmc_steps_evap);
    RUN_TEST(test_log_kmc_steps_evap_coord);

    RUN_TEST(test_write_xyz_suffix_iteration);
    RUN_TEST(test_write_xyz_suffix_time);

    // errors depend on implementation
    // RUN_TEST(test_fstring_to_buffer_returns_null_on_format_error);

    RUN_TEST(test_input_logging_basic);
    // RUN_TEST(test_write_xyz_file_creates_file_and_content);
    // RUN_TEST(test_write_logs_increments_framenums);
    // RUN_TEST(test_output_if_passed_checkpoint_triggers_write_logs);

    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
