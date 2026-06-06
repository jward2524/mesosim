#include "ErrorM.h"
#include "FileIO.h"
#include "Input.h"
#include "TUtils.h"
#include "Utils.h"
#include "unity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static struct SimulationState *ss = NULL;
static struct SimulationEnv *se = NULL;
static struct LoggingState *ls = NULL;

char *mock_name = "test/mock_input.in";
static FILE *mock_input_file = NULL;

void setUp(void)
{
    initialize_states(&ss, &se, &ls);
}

void tearDown(void)
{
    if (mock_input_file) {
        fclose(mock_input_file);
    };
    if (ss != NULL || se != NULL || ls != NULL) {
        clean_and_error(0);
    }
}

/* ================= Helpers ================= */

static FILE *open_mem(const char *text)
{
    FILE *tmpf = fopen(mock_name, "wb+");
    fopen_error(mock_name, tmpf);
    if (!tmpf) {
        perror("Failed to create temporary file");
        return NULL;
    }
    fputs(text, tmpf);
    rewind(tmpf);
    return tmpf;
}

static void assert_default_filename_time_deviation(char *time_str)
{
    time_t now = time(NULL);
    char *endptr;
    long file_time = strtol(time_str, &endptr, 10);
    TEST_ASSERT_TRUE_MESSAGE(*endptr == '\0', "Time part should be fully numeric");
    TEST_ASSERT_TRUE_MESSAGE(file_time > 0, "Extracted time from filename should be positive");

    long deviation = file_time - (long)now;
    if (deviation < 0) {
        deviation = -deviation;
    }
    long max_deviation = 1;
    char deviation_msg[64];
    snprintf(deviation_msg, sizeof(deviation_msg),
             "CSV filename time deviation should be <= %ld second(s)", max_deviation);
    TEST_ASSERT_TRUE_MESSAGE(deviation <= max_deviation, deviation_msg);
}

/* ================= Tests ================= */

/* === End to end tests === */

void test_parse_input_two_atomtypes_one_shell_success(void)
{
    const char *input = "atomtype Ag Au\n"
                        "composition 0.75 0.25\n"
                        "dissolution true false\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);
    finalize_atom_dependent(&ctx, &inputs);
    finalize_nne(&ctx, &inputs);

    TEST_ASSERT_EQUAL_INT_MESSAGE(2, inputs.atom_names_cnt, "Number of atom types");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, inputs.num_nn_levels, "Number of nn levels");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, inputs.substrate_composition[0],
                                     "Composition for atom type 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, inputs.substrate_composition[1],
                                     "Composition for atom type 1");

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, inputs.is_soluble[0], "Solubility for atom type 0");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, inputs.is_soluble[1], "Solubility for atom type 1");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.1, inputs.nn_energy[0], "1nne energy for atom type 0-0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.2, inputs.nn_energy[1], "1nne energy for atom type 0-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.3, inputs.nn_energy[2], "1nne energy for atom type 1-1");
}

void test_parse_input_three_atomtypes_two_shells_success(void)
{
    const char *input = "atomtype A B C\n"
                        "nnlevels 2\n"
                        "1nne 1 2 3 4 5 6\n"
                        "2nne 6 5 4 3 2 1\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);
    finalize_atom_dependent(&ctx, &inputs);
    finalize_nne(&ctx, &inputs);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, inputs.num_elements, "Number of elements");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, inputs.num_nn_levels, "Number of nn levels");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1, inputs.nn_energy[0],
                                     "1st shell nne energy for atom type 0-0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(2, inputs.nn_energy[1],
                                     "1st shell nne energy for atom type 0-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(3, inputs.nn_energy[2],
                                     "1st shell nne energy for atom type 0-2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(4, inputs.nn_energy[3],
                                     "1st shell nne energy for atom type 1-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(5, inputs.nn_energy[4],
                                     "1st shell nne energy for atom type 1-2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(6, inputs.nn_energy[5],
                                     "1st shell nne energy for atom type 2-2");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(6, inputs.nn_energy[6],
                                     "2nd shell nne energy for atom type 0-0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(5, inputs.nn_energy[7],
                                     "2nd shell nne energy for atom type 0-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(4, inputs.nn_energy[8],
                                     "2nd shell nne energy for atom type 0-2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(3, inputs.nn_energy[9],
                                     "2nd shell nne energy for atom type 1-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(2, inputs.nn_energy[10],
                                     "2nd shell nne energy for atom type 1-2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1, inputs.nn_energy[11],
                                     "2nd shell nne energy for atom type 2-2");
}

void test_parse_input_cluster_nns_file_success(void)
{
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_file("test/cluster_nns.in");
    TEST_ASSERT_NOT_NULL_MESSAGE(mock_input_file, "Failed to open test/cluster_nns.in");

    parse_input_file(mock_input_file, &ctx, &inputs, ls);
    finalize_config(&ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(128, inputs.system_size_x, "System size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, inputs.system_size_y, "System size y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, inputs.system_size_z, "System size z");

    TEST_ASSERT_EQUAL_INT_MESSAGE(GEOMETRY_CLUSTER, inputs.geometry, "Geometry type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(32, inputs.geometry_param, "Cluster radius");

    TEST_ASSERT_EQUAL_INT_MESSAGE(2, inputs.atom_names_cnt, "Number of atom types");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Ag", inputs.atom_names[0], "Atom type name at index 0");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Au", inputs.atom_names[1], "Atom type name at index 1");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, inputs.substrate_composition[0],
                                     "Composition for atom type 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, inputs.substrate_composition[1],
                                     "Composition for atom type 1");
    TEST_ASSERT_TRUE_MESSAGE(inputs.is_soluble[0], "Solubility for atom type 0");
    TEST_ASSERT_FALSE_MESSAGE(inputs.is_soluble[1], "Solubility for atom type 1");

    TEST_ASSERT_EQUAL_INT_MESSAGE(FCC, inputs.lattice_type, "Lattice type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, inputs.num_nn_levels, "Number of nn levels");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, inputs.nn_energy[0], "1st shell nne energy at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, inputs.nn_energy[1], "1st shell nne energy at index 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, inputs.nn_energy[2], "1st shell nne energy at index 2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, inputs.nn_energy[3], "2nd shell nne energy at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, inputs.nn_energy[4], "2nd shell nne energy at index 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, inputs.nn_energy[5], "2nd shell nne energy at index 2");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.9, inputs.initial_overpotential, "Initial overpotential");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.03, inputs.overpotential_ramp_rate,
                                     "Overpotential ramp rate");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.2, inputs.max_overpotential, "Maximum overpotential");

    TEST_ASSERT_EQUAL_INT_MESSAGE(2, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_CSV, ls->out_formats[0].type,
                                  "First output type should be CSV");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_XYZ, ls->out_formats[1].type,
                                  "Second output type should be XYZ");

    // Output CSV command
    OutputFormat *csv_format = &ls->out_formats[0];
    TEST_ASSERT_EQUAL_STRING_MESSAGE("test/output/cluster.csv", csv_format->csv.filename,
                                     "CSV filename");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_ITERATION, csv_format->csv.schedule.mode,
                                  "CSV schedule mode");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200.0, csv_format->csv.schedule.interval,
                                     "CSV interval value");
    // Check next_log_checkpoint for interval iteration mode
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200.0, csv_format->csv.schedule.next_checkpoint,
                                     "CSV next_csv_checkpoint should match first checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, csv_format->csv.field_count, "CSV field count");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("iter", csv_format->csv.field_names[0], "CSV field 0");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("time", csv_format->csv.field_names[1], "CSV field 1");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("energy", csv_format->csv.field_names[2], "CSV field 2");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("temperature", csv_format->csv.field_names[3], "CSV field 3");

    // Output XYZ command
    OutputFormat *xyz_format = &ls->out_formats[1];
    TEST_ASSERT_TRUE_MESSAGE(xyz_format->xyz.stripped, "XYZ stripped flag");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("test/output/cluster", xyz_format->xyz.prefix,
                                     "XYZ filename prefix");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_ITERATION, xyz_format->xyz.schedule.mode,
                                  "XYZ schedule mode");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(500, xyz_format->xyz.schedule.interval, "XYZ interval value");
    // Check next_log_checkpoint for interval time mode
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(500, xyz_format->xyz.schedule.next_checkpoint,
                                     "XYZ next_xyz_checkpoint should match first checkpoint");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293.0, inputs.temperature, "Temperature value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(12345U, inputs.rand_seed, "Random seed value");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FLAVOR_KMC, inputs.flavor, "Simulation flavor");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, inputs.sim_end_type,
                                  "Simulation end type");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(2000U, (unsigned int)inputs.final_iteration,
                                   "Final iteration value");
}

void test_parse_input_multi_command_unknown_mid_file_fails(void)
{
    const char *input = "systemsize 16 16 16\n"
                        "geometry cluster 4\n"
                        "atomtype A B\n"
                        "composition 0.5 0.5\n"
                        "bogus 123\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_UNKNOWN_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected unknown command failure in multi-command file");
    });
}

/* === Dependency checks === */

void test_parse_input_multi_command_finalize_composition_sum_fails(void)
{
    const char *input = "systemsize 16 16 16\n"
                        "geometry cluster 4\n"
                        "atomtype A B\n"
                        "composition 0.8 0.8\n"
                        "dissolution true false\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_DEPENDENCY, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_atom_dependent(&ctx, &inputs);
        finalize_nne(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected composition sum failure at finalize for multi-command file");
    });
}

void test_parse_input_missing_atomtype_fails(void)
{
    const char *input = "composition 1.0\n"
                        "nnlevels 1\n"
                        "1nne 0.1\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_DEPENDENCY, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_atom_dependent(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected failure");
    });
}

void test_parse_input_missing_nne_shell_fails(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 2\n"
                        "1nne 0.1 0.2 0.3\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_nne(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected failure");
    });
}

/* === Individual command checks === */

void test_parse_input_unknown_command_fails(void)
{
    const char *input = "atomtype A B\n"
                        "foobar 123\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_UNKNOWN_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected an unknown command failure");
    });
}

void test_nnlevels_wrong_argcount_fails(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 1 2\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected nnlevels argcount failure");
    });
}

void test_composition_non_numeric_fails(void)
{
    const char *input = "atomtype A B\n"
                        "composition 0.5 foo\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected non-numeric composition failure");
    });
}

void test_dissolution_invalid_boolean_fails(void)
{
    const char *input = "atomtype A B\n"
                        "dissolution true maybe\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected invalid boolean failure");
    });
}

void test_nne_level_exceeds_nnlevels_fails(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "2nne 0.1 0.2 0.3\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_nne(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected nne level overflow failure");
    });
}

void test_nne_duplicate_definition_fails(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n"
                        "1nne 0.3 0.2 0.1\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_DUPLICATE_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_nne(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected duplicate nne failure");
    });
}

void test_nne_too_few_values_fails(void)
{
    const char *input = "atomtype A B C\n"
                        "nnlevels 1\n"
                        "1nne 1 2 3\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_nne(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected nne value count failure");
    });
}

void test_composition_count_mismatch_fails(void)
{
    const char *input = "atomtype A B C\n"
                        "composition 0.5 0.5\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_atom_dependent(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected composition count mismatch");
    });
}

void test_systemsize_valid_values_success(void)
{
    const char *input = "systemsize 10 11 12\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(10, inputs.system_size_x, "System size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(11, inputs.system_size_y, "System size y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12, inputs.system_size_z, "System size z");
}

void test_systemsize_non_numeric_fails(void)
{
    const char *input = "systemsize 10 a 12\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected systemsize parse failure");
    });
}

void test_temp_valid_value_success(void)
{
    const char *input = "temp 293.5\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293.5, inputs.temperature, "Temperature value");
}

void test_temp_non_numeric_fails(void)
{
    const char *input = "temp abc\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected temp parse failure");
    });
}

void test_seed_default_keyword_success(void)
{
    const char *input = "seed default\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(DEFAULT_SEED, inputs.rand_seed, "Default random seed");
}

void test_seed_non_numeric_fails(void)
{
    const char *input = "seed notanumber\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected seed parse failure");
    });
}

void test_potential_sweep_values_success(void)
{
    const char *input = "potential 0.9 0.03 1.2\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.9, inputs.initial_overpotential, "Initial overpotential");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.03, inputs.overpotential_ramp_rate,
                                     "Overpotential ramp rate");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.2, inputs.max_overpotential, "Maximum overpotential");
}

void test_potential_bad_argcount_fails(void)
{
    const char *input = "potential 0.9 0.03\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected potential argcount failure");
    });
}

void test_struct_fcc_success(void)
{
    const char *input = "struct FCC\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(FCC, inputs.lattice_type, "Lattice type");
}

void test_struct_invalid_value_fails(void)
{
    const char *input = "struct HCP\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected struct value failure");
    });
}

void test_output_csv_interval_iteration_success(void)
{
    const char *input = "output csv test/output/cluster.csv interval iteration 200 fields iter "
                        "time energy temperature\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_CSV, ls->out_formats[0].type,
                                  "First output type should be CSV");

    OutputFormat *format = &ls->out_formats[0];
    TEST_ASSERT_EQUAL_STRING_MESSAGE("test/output/cluster.csv", format->csv.filename,
                                     "CSV filename");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_ITERATION, format->csv.schedule.mode,
                                  "CSV schedule iteration interval mode");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200.0, format->csv.schedule.interval, "CSV interval value");
    // Check next_log_checkpoint for interval iteration mode
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200.0, format->csv.schedule.next_checkpoint,
                                     "CSV next_csv_checkpoint should match first checkpoint");
    TEST_ASSERT_NULL_MESSAGE(format->csv.schedule.list, "CSV list value in interval mode");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, format->csv.field_count, "CSV field count");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("iter", format->csv.field_names[0], "CSV field 0");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("time", format->csv.field_names[1], "CSV field 1");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("energy", format->csv.field_names[2], "CSV field 2");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("temperature", format->csv.field_names[3], "CSV field 3");
}

void test_output_csv_list_time_success(void)
{
    const char *input =
        "output csv test/output/cluster.csv list time 0.1 0.2 0.3 fields iter time energy\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_CSV, ls->out_formats[0].type,
                                  "First output type should be CSV");

    OutputFormat *format = &ls->out_formats[0];
    TEST_ASSERT_EQUAL_STRING_MESSAGE("test/output/cluster.csv", format->csv.filename,
                                     "CSV filename");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_LIST_TIME, format->csv.schedule.mode,
                                  "CSV list time schedule mode");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, format->csv.schedule.list_len, "CSV schedule list length");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0, format->csv.schedule.interval,
                                     "CSV interval value in list mode");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.1, format->csv.schedule.list[0],
                                     "CSV schedule list value 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.2, format->csv.schedule.list[1],
                                     "CSV schedule list value 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.3, format->csv.schedule.list[2],
                                     "CSV schedule list value 2");
    // Check next_log_checkpoint for list time mode
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.1, format->csv.schedule.next_checkpoint,
                                     "CSV next_checkpoint should match first checkpoint in list");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, format->csv.field_count, "CSV field count");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("iter", format->csv.field_names[0], "CSV field 0");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("time", format->csv.field_names[1], "CSV field 1");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("energy", format->csv.field_names[2], "CSV field 2");
}

void test_output_xyz_interval_time_success(void)
{
    const char *input = "output xyz test/output/cluster interval time 0.5\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_XYZ, ls->out_formats[0].type,
                                  "First output type should be XYZ");

    OutputFormat *format = &ls->out_formats[0];
    TEST_ASSERT_EQUAL_STRING_MESSAGE("test/output/cluster", format->xyz.prefix,
                                     "XYZ filename prefix");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_TIME, format->xyz.schedule.mode,
                                  "XYZ schedule interval time mode");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.5, format->xyz.schedule.interval, "XYZ interval value");
    // Check next_log_checkpoint for interval time mode
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.5, format->xyz.schedule.next_checkpoint,
                                     "XYZ next_xyz_checkpoint should match first checkpoint");
    TEST_ASSERT_FALSE_MESSAGE(format->xyz.stripped, "XYZ stripped flag should be false by default");
}

void test_output_xyz_stripped_success(void)
{
    const char *input = "output xyz stripped test/output/cluster interval time 0.2\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_XYZ, ls->out_formats[0].type,
                                  "First output type should be XYZ");

    OutputFormat *format = &ls->out_formats[0];
    TEST_ASSERT_TRUE_MESSAGE(format->xyz.stripped,
                             "XYZ stripped flag should be true when 'stripped' is present");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("test/output/cluster", format->xyz.prefix,
                                     "XYZ filename prefix");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_SCHEDULE_INTERVAL_TIME, format->xyz.schedule.mode,
                                  "XYZ schedule interval time mode");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.2, format->xyz.schedule.interval, "XYZ interval value");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.2, format->xyz.schedule.next_checkpoint,
                                     "XYZ next_xyz_checkpoint should match first checkpoint");
}

void test_output_iter_default_filename_success(void)
{
    const char *input = "output steps\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_STEPS_CSV, ls->out_formats[0].type,
                                  "First output type should be a steps csv");

    OutputFormat *format = &ls->out_formats[0];
    TEST_ASSERT_TRUE_MESSAGE(strlen(format->steps.filename) > 0,
                             "Default steps filename should be set");
    // Should end with _steps.csv
    size_t len = strlen(format->steps.filename);
    TEST_ASSERT_TRUE_MESSAGE(len > 10 &&
                                 strcmp(format->steps.filename + len - 10, "_steps.csv") == 0,
                             "Default steps filename should end with _steps.csv");
}

void test_output_steps_default_filename_coord_success(void)
{
    const char *input = "output steps coord\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_STEPS_CSV, ls->out_formats[0].type,
                                  "First output type should be a steps csv");

    OutputFormat *format = &ls->out_formats[0];
    TEST_ASSERT_TRUE_MESSAGE(strlen(format->steps.filename) > 0,
                             "Default steps filename should be set");
    TEST_ASSERT_TRUE_MESSAGE(format->steps.with_coordination, "Steps CSV coordination");
    // Should end with _steps.csv
    size_t len = strlen(format->steps.filename);
    TEST_ASSERT_TRUE_MESSAGE(len > 10 &&
                                 strcmp(format->steps.filename + len - 10, "_steps.csv") == 0,
                             "Default steps filename should end with _steps.csv");
}

void test_output_steps_with_filename_success(void)
{
    const char *input = "output steps my_steps.csv\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_STEPS_CSV, ls->out_formats[0].type,
                                  "First output type should be a steps csv");

    OutputFormat *format = &ls->out_formats[0];
    TEST_ASSERT_EQUAL_STRING_MESSAGE("my_steps.csv", format->steps.filename,
                                     "Steps filename should match provided filename");
}

void test_output_steps_with_filename_coord_success(void)
{
    const char *input = "output steps my_steps.csv coord\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_STEPS_CSV, ls->out_formats[0].type,
                                  "First output type should be a steps csv");

    OutputFormat *format = &ls->out_formats[0];
    TEST_ASSERT_EQUAL_STRING_MESSAGE("my_steps.csv", format->steps.filename,
                                     "Steps filename should match provided filename");
    TEST_ASSERT_TRUE_MESSAGE(format->steps.with_coordination, "Steps csv coord flag");
}

void test_output_steps_with_filename_and_extra_args_fails(void)
{
    const char *input = "output steps my_steps.csv extra\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected failure for 'output iter' with filename and extra arguments");
    });
}

void test_output_invalid_mode_fails(void)
{
    const char *input =
        "output csv test/output/cluster.csv cadence iteration 200 fields iter time\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected output mode failure");
    });
}

void test_output_missing_fields_keyword_fails(void)
{
    const char *input =
        "output csv test/output/cluster.csv interval iteration 200 iter time energy\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected missing fields keyword failure");
    });
}

void test_output_fields_empty_fails(void)
{
    const char *input = "output csv test/output/cluster.csv interval iteration 200 fields\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected empty fields failure");
    });
}

void test_output_xyz_with_fields_fails(void)
{
    const char *input =
        "output xyz test/output/cluster.xyz interval iteration 200 fields iter time\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected failure for output xyz with fields");
    });
}

void test_output_csv_default_filename_time_success(void)
{
    const char *input = "output csv interval iteration 100 fields iter time\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_CSV, ls->out_formats[0].type,
                                  "First output type should be CSV");

    OutputFormat *format = &ls->out_formats[0];

    // Extract the time from the filename, which should be of the form "[time].csv"
    TEST_ASSERT_NOT_NULL_MESSAGE(format->csv.filename, "CSV filename should not be NULL");
    const char *dot = strrchr(format->csv.filename, '.');
    TEST_ASSERT_NOT_NULL_MESSAGE(dot, "CSV filename should contain a dot");
    TEST_ASSERT_EQUAL_STRING_MESSAGE(".csv", dot, "CSV filename should end with .csv");

    size_t len = (size_t)(dot - format->csv.filename);
    char time_part[32];
    TEST_ASSERT_TRUE_MESSAGE(len < sizeof(time_part), "Time part too long");
    strncpy(time_part, format->csv.filename, len);
    time_part[len] = '\0';

    assert_default_filename_time_deviation(time_part);
}

void test_output_xyz_default_prefix_time_success(void)
{
    const char *input = "output xyz interval iteration 100\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);
    time_t now = time(NULL);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls->out_formats_cnt, "Number of output formats");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls->out_formats, "Output formats array should not be null");
    TEST_ASSERT_EQUAL_INT_MESSAGE(OUTPUT_FORMAT_XYZ, ls->out_formats[0].type,
                                  "First output type should be XYZ");

    OutputFormat *format = &ls->out_formats[0];

    // Extract the time from the prefix, which should be of the form "[time].xyz"
    TEST_ASSERT_NOT_NULL_MESSAGE(format->xyz.prefix, "XYZ prefix should not be NULL");
    const char *dot = strrchr(format->xyz.prefix, '.');
    TEST_ASSERT_NOT_NULL_MESSAGE(dot, "XYZ prefix should contain a dot");
    size_t len = (size_t)(dot - format->xyz.prefix);
    char time_part[32];
    TEST_ASSERT_TRUE_MESSAGE(len < sizeof(time_part), "Time part too long");
    strncpy(time_part, format->xyz.prefix, len);
    time_part[len] = '\0';
    char *endptr;
    long file_time = strtol(time_part, &endptr, 10);
    TEST_ASSERT_TRUE_MESSAGE(*endptr == '\0', "XYZ prefix time part should be fully numeric");
    TEST_ASSERT_TRUE_MESSAGE(file_time > 0, "Extracted time from prefix should be positive");

    long deviation = file_time - (long)now;
    if (deviation < 0) {
        deviation = -deviation;
    }
    long max_deviation = 1;
    char deviation_msg[64];
    snprintf(deviation_msg, sizeof(deviation_msg),
             "XYZ prefix time deviation should be <= %ld second(s)", max_deviation);
    TEST_ASSERT_TRUE_MESSAGE(deviation <= max_deviation, deviation_msg);
    TEST_ASSERT_FALSE_MESSAGE(format->xyz.stripped, "XYZ stripped flag should be false by default");
}

void test_output_unrecognized_field_fails(void)
{
    const char *input =
        "output csv test/output/cluster.csv interval iteration 200 fields iter time notafield\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected failure for unrecognized field");
    });
}

void test_output_unsupported_field_for_flavor_fails(void)
{
    const char *input =
        "flavor KMC\n"
        "output csv test/output/cluster.csv interval iteration 200 fields iter time notafield\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_csv_fields(ls, inputs.flavor);
        TEST_FAIL_MESSAGE("Expected failure for unrecognized field");
    });
}

void test_output_non_numeric_interval_fails(void)
{
    const char *input =
        "output csv test/output/cluster.csv interval iteration notanumber fields iter time\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected failure for non-numeric interval");
    });
}

void test_checkpoint_interval_only_success(void)
{
    const char *input = "checkpoint 500\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(500, ls->checkpoint.interval, "Checkpoint schedule interval");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(500, ls->checkpoint.next_checkpoint,
                                     "Checkpoint next_checkpoint should match interval");

    // check if default filename is correct
    TEST_ASSERT_TRUE_MESSAGE(strlen(ls->checkpoint.filename) > 0,
                             "Default checkpoint filename should be set");

    const char *dot = strrchr(ls->checkpoint.filename, '.');
    TEST_ASSERT_NOT_NULL_MESSAGE(dot, "Checkpoint filename should contain a dot");
    TEST_ASSERT_EQUAL_STRING_MESSAGE(".bin", dot, "Checkpoint filename should end with .bin");

    const char *under = strrchr(ls->checkpoint.filename, '_');
    TEST_ASSERT_NOT_NULL_MESSAGE(under, "Checkpoint filename should contain an underscore");
    char check_prefix[32];
    size_t prefix_len = (size_t)(under - ls->checkpoint.filename);
    TEST_ASSERT_TRUE_MESSAGE(prefix_len < sizeof(check_prefix), "Checkpoint prefix too long");
    strncpy(check_prefix, ls->checkpoint.filename, prefix_len);
    check_prefix[prefix_len] = '\0';
    TEST_ASSERT_EQUAL_STRING_MESSAGE("checkpoint", check_prefix,
                                     "Checkpoint filename should start with 'checkpoint'");

    // Extract the time from the filename, which should be of the form "checkpoint_[time].bin"
    const char *time_start = under + 1;
    char time_part[32];
    size_t time_len = (size_t)(dot - time_start);
    TEST_ASSERT_TRUE_MESSAGE(time_len < sizeof(time_part), "Checkpoint time part too long");
    strncpy(time_part, time_start, time_len);
    time_part[time_len] = '\0';
    assert_default_filename_time_deviation(time_part);
}

void test_geometry_cluster_success(void)
{
    const char *input = "geometry cluster 32\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(GEOMETRY_CLUSTER, inputs.geometry, "Geometry type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(32, inputs.geometry_param, "Cluster radius");
}

void test_geometry_invalid_type_fails(void)
{
    const char *input = "geometry cube 4\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected geometry type failure");
    });
}

void test_run_iteration_mode_success(void)
{
    const char *input = "run iteration 2000\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, inputs.sim_end_type,
                                  "Simulation end type");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(2000U, (unsigned int)inputs.final_iteration,
                                   "Final iteration value");
}

void test_run_unknown_mode_fails(void)
{
    const char *input = "run steps 100\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected run mode failure");
    });
}

void test_flavor_kmc_success(void)
{
    const char *input = "flavor KMC\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(FLAVOR_KMC, inputs.flavor, "Simulation flavor");
}

void test_flavor_invalid_value_fails(void)
{
    const char *input = "flavor BAD\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected flavor value failure");
    });
}

void test_systemsize_wrong_argcount_fails(void)
{
    const char *input = "systemsize 10 11\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected systemsize argcount failure");
    });
}

void test_temp_extra_arg_fails(void)
{
    const char *input = "temp 300 301\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected temp argcount failure");
    });
}

void test_seed_missing_arg_fails(void)
{
    const char *input = "seed\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected seed missing arg failure");
    });
}

void test_potential_non_numeric_sweep_fails(void)
{
    const char *input = "potential 0.9 foo 1.2\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected potential numeric parse failure");
    });
}

void test_geometry_cluster_non_numeric_fails(void)
{
    const char *input = "geometry cluster radius\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected geometry cluster parse failure");
    });
}

void test_geometry_file_missing_name_fails(void)
{
    const char *input = "geometry file\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected geometry file missing name failure");
    });
}

void test_run_iteration_non_numeric_fails(void)
{
    const char *input = "run iteration ten\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected run iteration parse failure");
    });
}

void test_run_bad_argcount_fails(void)
{
    const char *input = "run time\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected run argcount failure");
    });
}

void test_nne_non_numeric_value_fails(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne 0.1 x 0.3\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected nne numeric parse failure");
    });
}

void test_composition_sum_not_one_fails(void)
{
    const char *input = "atomtype A B\n"
                        "composition 0.6 0.5\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_DEPENDENCY, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_atom_dependent(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected composition sum finalize failure");
    });
}

void test_dissolution_count_mismatch_additional_fails(void)
{
    const char *input = "atomtype A B C\n"
                        "dissolution true false\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_atom_dependent(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected dissolution count mismatch finalize failure");
    });
}

void test_atomtype_three_types_success(void)
{
    const char *input = "atomtype Ag Au Cu\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, inputs.atom_names_cnt, "Number of atom names");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, inputs.num_elements, "Number of elements");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Ag", inputs.atom_names[0], "Atom type name at index 0");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Au", inputs.atom_names[1], "Atom type name at index 1");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Cu", inputs.atom_names[2], "Atom type name at index 2");
}

void test_atomtype_missing_args_fails(void)
{
    const char *input = "atomtype\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected atomtype missing-args failure");
    });
}

void test_composition_two_types_finalize_success(void)
{
    const char *input = "atomtype A B\n"
                        "composition 0.25 0.75\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);
    finalize_atom_dependent(&ctx, &inputs);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, inputs.substrate_composition[0],
                                     "Composition at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, inputs.substrate_composition[1],
                                     "Composition at index 1");
}

void test_composition_missing_values_fails(void)
{
    const char *input = "composition\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected composition missing-values failure");
    });
}

void test_dissolution_all_true_sets_flag_success(void)
{
    const char *input = "atomtype A B C\n"
                        "dissolution true true true\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);
    finalize_atom_dependent(&ctx, &inputs);

    TEST_ASSERT_TRUE_MESSAGE(inputs.is_soluble[0], "Solubility at index 0");
    TEST_ASSERT_TRUE_MESSAGE(inputs.is_soluble[1], "Solubility at index 1");
    TEST_ASSERT_TRUE_MESSAGE(inputs.is_soluble[2], "Solubility at index 2");
}

void test_dissolution_missing_values_fails(void)
{
    const char *input = "dissolution\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected dissolution missing-values failure");
    });
}

void test_nnlevels_single_value_success(void)
{
    const char *input = "nnlevels 3\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, inputs.num_nn_levels, "NN levels value");
}

void test_nnlevels_non_numeric_fails(void)
{
    const char *input = "nnlevels abc\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected nnlevels non-numeric failure");
    });
}

void test_nne_missing_values_fails(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_COUNT_MISMATCH, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        TEST_FAIL_MESSAGE("Expected nne missing-values failure");
    });
}

/* === Finalization checks === */

void test_finalize_nne_missing_nnlevels_fails(void)
{
    const char *input = "atomtype A B\n"
                        "1nne 0.1 0.2 0.3\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_DEPENDENCY, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_nne(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected finalize_nne without nnlevels failure");
    });
}

void test_finalize_atom_dependent_direct_success(void)
{
    const char *input = "atomtype A B\n"
                        "composition 0.4 0.6\n"
                        "dissolution false true\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);
    finalize_atom_dependent(&ctx, &inputs);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.4, inputs.substrate_composition[0],
                                     "Finalized composition at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.6, inputs.substrate_composition[1],
                                     "Finalized composition at index 1");
    TEST_ASSERT_FALSE_MESSAGE(inputs.is_soluble[0], "Finalized solubility at index 0");
    TEST_ASSERT_TRUE_MESSAGE(inputs.is_soluble[1], "Finalized solubility at index 1");
}

void test_finalize_atom_dependent_missing_atomtype_fails(void)
{
    const char *input = "composition 1.0\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_DEPENDENCY, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_atom_dependent(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected finalize_atom_dependent missing atomtype failure");
    });
}

void test_finalize_nne_direct_single_shell_success(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne 0.11 0.22 0.33\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);
    finalize_nne(&ctx, &inputs);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, inputs.num_nn_types, "Total NN energy types");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.11, inputs.nn_energy[0], "NN energy at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.22, inputs.nn_energy[1], "NN energy at index 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.33, inputs.nn_energy[2], "NN energy at index 2");
}

void test_finalize_nne_direct_duplicate_shell_fails(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n"
                        "1nne 0.3 0.2 0.1\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_DUPLICATE_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_nne(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected finalize_nne duplicate shell failure");
    });
}

void test_finalize_nne_direct_level_exceeds_nnlevels_fails(void)
{
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "2nne 0.1 0.2 0.3\n";
    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_nne(&ctx, &inputs);
        TEST_FAIL_MESSAGE("Expected finalize_nne level overflow failure");
    });
}

/* === Required commands checks === */
void test_required_commands_all_present_success(void)
{
    const char *input = "systemsize 8 8 8\n"
                        "struct FCC\n"
                        "geometry cluster 4\n"
                        "atomtype Ag Au\n"
                        "composition 0.5 0.5\n"
                        "dissolution true false\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n"
                        "run iteration 10\n"
                        "flavor KMC\n"
                        "temp 300\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    parse_input_file(mock_input_file, &ctx, &inputs, ls);
    finalize_config(&ctx, &inputs, ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(8, inputs.system_size_x, "System size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FCC, inputs.lattice_type, "Lattice type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(GEOMETRY_CLUSTER, inputs.geometry, "Geometry type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, inputs.atom_names_cnt, "Atom type count");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.5, inputs.substrate_composition[0], "Composition 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.5, inputs.substrate_composition[1], "Composition 1");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, inputs.is_soluble[0], "Solubility 0");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, inputs.is_soluble[1], "Solubility 1");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, inputs.num_nn_levels, "NN levels");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.1, inputs.nn_energy[0], "NN energy 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.2, inputs.nn_energy[1], "NN energy 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.3, inputs.nn_energy[2], "NN energy 2");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, inputs.sim_end_type, "Sim end type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(10, inputs.final_iteration, "Final iteration");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FLAVOR_KMC, inputs.flavor, "Flavor");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(300.0, inputs.temperature, "Temperature");
}

void test_missing_required_systemsize_fails(void)
{
    const char *input =
        /* "systemsize 8 8 8\n" intentionally omitted */
        "struct FCC\n"
        "geometry cluster 4\n"
        "atomtype Ag Au\n"
        "composition 0.5 0.5\n"
        "dissolution true false\n"
        "nnlevels 1\n"
        "1nne 0.1 0.2 0.3\n"
        "run iteration 10\n"
        "flavor KMC\n"
        "temp 300\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_config(&ctx, &inputs, ls);
        check_required_inputs(&ctx, inputs.flavor);
        TEST_FAIL_MESSAGE("Expected failure due to missing required command: systemsize");
    });
}

void test_missing_required_struct_fails(void)
{
    const char *input = "systemsize 8 8 8\n"
                        /* "struct FCC\n" intentionally omitted */
                        "geometry cluster 4\n"
                        "atomtype Ag Au\n"
                        "composition 0.5 0.5\n"
                        "dissolution true false\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n"
                        "run iteration 10\n"
                        "flavor KMC\n"
                        "temp 300\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_config(&ctx, &inputs, ls);
        check_required_inputs(&ctx, inputs.flavor);
        TEST_FAIL_MESSAGE("Expected failure due to missing required command: struct");
    });
}

void test_missing_multiple_required_commands_fails(void)
{
    const char *input =
        /* "systemsize 8 8 8\n" intentionally omitted */
        /* "struct FCC\n" intentionally omitted */
        "geometry cluster 4\n"
        "atomtype Ag Au\n"
        "composition 0.5 0.5\n"
        "dissolution true false\n"
        "nnlevels 1\n"
        "1nne 0.1 0.2 0.3\n"
        "run iteration 10\n"
        "flavor KMC\n"
        "temp 300\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_config(&ctx, &inputs, ls);
        check_required_inputs(&ctx, inputs.flavor);
        TEST_FAIL_MESSAGE(
            "Expected failure due to multiple missing required commands: systemsize, struct");
    });
}

void test_missing_all_required_commands_fails(void)
{
    const char *input = "# No required commands present\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        // finalize_config(&ctx, &inputs, ls);
        check_required_inputs(&ctx, inputs.flavor);
        TEST_FAIL_MESSAGE("Expected failure due to all required commands missing");
    });
}

void test_required_commands_invalid_args_fail_for_argument(void)
{
    const char *input = "systemsize 8 X 8\n" // 'X' is invalid
                        "struct FCC\n"
                        "geometry cluster 4\n"
                        "atomtype Ag Au\n"
                        "composition 0.5 0.5\n"
                        "dissolution true false\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n"
                        "run iteration 10\n"
                        "flavor KMC\n"
                        "temp 300\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_INVALID_ARG, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_config(&ctx, &inputs, ls);
        check_required_inputs(&ctx, inputs.flavor);
        TEST_FAIL_MESSAGE("Expected failure due to invalid argument in required command");
    });
}

void test_required_commands_commented_out_treated_as_missing(void)
{
    const char *input = "# systemsize 8 8 8\n" // commented out
                        "struct FCC\n"
                        "geometry cluster 4\n"
                        "atomtype Ag Au\n"
                        "composition 0.5 0.5\n"
                        "dissolution true false\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n"
                        "run iteration 10\n"
                        "flavor KMC\n"
                        "temp 300\n";

    ParseContext ctx = {0};
    struct SimulationConfig inputs = {0};
    mock_input_file = open_mem(input);

    EXPECT_EXIT(INPUT_ERR_MISSING_CMD, {
        parse_input_file(mock_input_file, &ctx, &inputs, ls);
        finalize_config(&ctx, &inputs, ls);
        check_required_inputs(&ctx, inputs.flavor);
        TEST_FAIL_MESSAGE("Expected failure due to required command being commented out");
    });
}

/* ================= Runner ================= */

int main(void)
{
    UNITY_BEGIN();

    // Command requirements
    RUN_TEST(test_required_commands_invalid_args_fail_for_argument);
    RUN_TEST(test_required_commands_commented_out_treated_as_missing);

    // atomtype
    RUN_TEST(test_atomtype_three_types_success);
    RUN_TEST(test_atomtype_missing_args_fails);

    // composition
    RUN_TEST(test_composition_two_types_finalize_success);
    RUN_TEST(test_composition_non_numeric_fails);
    RUN_TEST(test_composition_count_mismatch_fails);
    RUN_TEST(test_composition_missing_values_fails);
    RUN_TEST(test_composition_sum_not_one_fails);

    // dissolution
    RUN_TEST(test_dissolution_all_true_sets_flag_success);
    RUN_TEST(test_dissolution_invalid_boolean_fails);
    RUN_TEST(test_dissolution_missing_values_fails);
    RUN_TEST(test_dissolution_count_mismatch_additional_fails);

    // nnlevels
    RUN_TEST(test_nnlevels_single_value_success);
    RUN_TEST(test_nnlevels_wrong_argcount_fails);
    RUN_TEST(test_nnlevels_non_numeric_fails);

    // nne
    RUN_TEST(test_nne_level_exceeds_nnlevels_fails);
    RUN_TEST(test_nne_duplicate_definition_fails);
    RUN_TEST(test_nne_too_few_values_fails);
    RUN_TEST(test_nne_non_numeric_value_fails);
    RUN_TEST(test_nne_missing_values_fails);

    // systemsize
    RUN_TEST(test_systemsize_valid_values_success);
    RUN_TEST(test_systemsize_non_numeric_fails);
    RUN_TEST(test_systemsize_wrong_argcount_fails);

    // temp
    RUN_TEST(test_temp_valid_value_success);
    RUN_TEST(test_temp_non_numeric_fails);
    RUN_TEST(test_temp_extra_arg_fails);

    // seed
    RUN_TEST(test_seed_default_keyword_success);
    RUN_TEST(test_seed_non_numeric_fails);
    RUN_TEST(test_seed_missing_arg_fails);

    // potential
    RUN_TEST(test_potential_sweep_values_success);
    RUN_TEST(test_potential_bad_argcount_fails);
    RUN_TEST(test_potential_non_numeric_sweep_fails);

    // struct
    RUN_TEST(test_struct_fcc_success);
    RUN_TEST(test_struct_invalid_value_fails);

    // output
    RUN_TEST(test_output_csv_interval_iteration_success);
    RUN_TEST(test_output_csv_list_time_success);
    RUN_TEST(test_output_xyz_interval_time_success);
    RUN_TEST(test_output_iter_default_filename_success);
    RUN_TEST(test_output_steps_with_filename_success);
    RUN_TEST(test_output_steps_with_filename_and_extra_args_fails);
    RUN_TEST(test_output_invalid_mode_fails);
    RUN_TEST(test_output_missing_fields_keyword_fails);
    RUN_TEST(test_output_fields_empty_fails);
    RUN_TEST(test_output_xyz_with_fields_fails);
    RUN_TEST(test_output_unrecognized_field_fails);
    RUN_TEST(test_output_non_numeric_interval_fails);
    RUN_TEST(test_output_csv_default_filename_time_success);
    RUN_TEST(test_output_xyz_default_prefix_time_success);
    RUN_TEST(test_output_xyz_stripped_success);

    // checkpoint
    RUN_TEST(test_checkpoint_interval_only_success);
    // RUN_TEST(test_checkpoint_filename_success);
    // RUN_TEST(test_checkpoint_filename_only_fails);
    // RUN_TEST(test_checkpoint_no_arguments_fails);

    // geometry
    RUN_TEST(test_geometry_cluster_success);
    RUN_TEST(test_geometry_invalid_type_fails);
    RUN_TEST(test_geometry_cluster_non_numeric_fails);
    RUN_TEST(test_geometry_file_missing_name_fails);

    // run
    RUN_TEST(test_run_iteration_mode_success);
    RUN_TEST(test_run_unknown_mode_fails);
    RUN_TEST(test_run_iteration_non_numeric_fails);
    RUN_TEST(test_run_bad_argcount_fails);

    // flavor
    RUN_TEST(test_flavor_kmc_success);
    RUN_TEST(test_flavor_invalid_value_fails);

    // Direct finalizer tests
    RUN_TEST(test_finalize_nne_missing_nnlevels_fails);
    RUN_TEST(test_finalize_atom_dependent_direct_success);
    RUN_TEST(test_finalize_atom_dependent_missing_atomtype_fails);
    RUN_TEST(test_finalize_nne_direct_single_shell_success);
    RUN_TEST(test_finalize_nne_direct_duplicate_shell_fails);
    RUN_TEST(test_finalize_nne_direct_level_exceeds_nnlevels_fails);

    // End-to-end parser baselines
    RUN_TEST(test_parse_input_two_atomtypes_one_shell_success);
    RUN_TEST(test_parse_input_three_atomtypes_two_shells_success);

    RUN_TEST(test_parse_input_cluster_nns_file_success);
    RUN_TEST(test_parse_input_multi_command_unknown_mid_file_fails);

    RUN_TEST(test_parse_input_multi_command_finalize_composition_sum_fails);
    RUN_TEST(test_parse_input_missing_atomtype_fails);
    RUN_TEST(test_parse_input_missing_nne_shell_fails);
    RUN_TEST(test_parse_input_unknown_command_fails);

    RUN_TEST(test_required_commands_all_present_success);
    RUN_TEST(test_missing_required_systemsize_fails);
    RUN_TEST(test_missing_required_struct_fails);
    RUN_TEST(test_missing_multiple_required_commands_fails);
    RUN_TEST(test_missing_all_required_commands_fails);

    UNITY_END();

    remove(mock_name);

    // return 0 else makefile throws error
    return 0;
}
