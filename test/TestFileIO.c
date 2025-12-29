#include "Atoms.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "State.h"
#include "Utils.h"
#include "unity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct SimulationState *ss;
struct SimulationEnv *se;
struct LoggingState *ls;
FILE *temp_log;
FILE *atom_file;
const char temp_name[] = "temp.log";

static void fopen_error(const char *filename, FILE *file)
{
    if (file == NULL) {
        printf("ERROR! Couldn't open output file %s\n", filename);
        fprintf(stderr, "Couldn't open file %s: %s\n", filename, strerror(errno));
        TEST_ASSERT_NOT_NULL_MESSAGE(file, "File not opened - check result file");
    }
}

// before and after each test (each RUN_TEST)
void setUp(void)
{
    ss = calloc(1, sizeof(struct SimulationState));
    se = calloc(1, sizeof(struct SimulationEnv));
    ls = calloc(1, sizeof(struct LoggingState));

    set_state(ss, se, ls);
    temp_log = fopen(temp_name, "w");
    fopen_error(temp_name, temp_log);
}

void tearDown(void)
{
    free(ss);
    free(se);
    free(ls);

    // fclose needs to be here in case a test fails
    fclose(temp_log);
    if (atom_file) {
        fclose(atom_file);
    }
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

    TEST_ASSERT_EQUAL_INT_MESSAGE(SIMULATION_TYPE_FROM_FILE, se->simulation_type, "system size x");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("test/sheet.xyz", se->atoms_filename, "system size x");
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
    char **arr = (char *[]){"Ag", "Au"};
    se->atom_names = arr;
    se->atom_names_cnt = 2;
    get_shifts(se);
    set_primitive_basis(se);
    initialize_simulation_box(se);
    initialize_zones(ss->zone_arr, se);
    initialize_simulation_variables(ss, se);

    process_xyz_file(temp_log, atom_file, ss, se, ls);

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

int main(void)
{
    UNITY_BEGIN();

    // tests to run
    RUN_TEST(test_parse_input_systemsize);
    RUN_TEST(test_parse_input_geometry_file);
    RUN_TEST(test_process_xyz_file);

    UNITY_END();

    if (temp_log) {
        int rc = remove(temp_name);
        if (rc)
            perror("Remove of test log file failed");
    }

    // return 0 else makefile throws error
    return 0;
}
