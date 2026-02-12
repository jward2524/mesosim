#include "unity.h"
#include "Input.h"
#include "FileIO.h"
#include "Utils.h"
#include "ErrorM.h"
#include "TUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static struct SimulationState *ss = NULL;
static struct SimulationEnv *se = NULL;
static struct LoggingState *ls = NULL;

/* ================= Helpers ================= */

static FILE *open_mem(const char *text) {
    FILE *tmpf = tmpfile();
    if (!tmpf) {
        perror("Failed to create temporary file");
        return NULL;
    }
    fputs(text, tmpf);
    rewind(tmpf);
    return tmpf;
}

void setUp(void)
{
    initialize_states(&ss, &se, &ls);
}

void tearDown(void)
{
    clean_and_error(0);
}

/* ================= Tests ================= */

void test_parse_input_two_atomtypes_one_shell_success(void) {
    const char *input = "atomtype Ag Au\n"
                        "composition 0.75 0.25\n"
                        "dissolution true false\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    finalize_config(&ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->atom_names_cnt, "Number of atom types");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->num_nn_levels, "Number of nn levels");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, se->substrate_composition[0], "Composition for atom type 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE( 0.25, se->substrate_composition[1], "Composition for atom type 1");

    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->is_soluble[0], "Solubility for atom type 0");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, se->is_soluble[1], "Solubility for atom type 1");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.1, se->nn_energy[0], "1nne energy for atom type 0-0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.2, se->nn_energy[1], "1nne energy for atom type 0-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.3, se->nn_energy[2], "1nne energy for atom type 1-1");
}

void test_parse_input_three_atomtypes_two_shells_success(void) {
    const char *input = "atomtype A B C\n"
                        "nnlevels 2\n"
                        "1nne 1 2 3 4 5 6\n"
                        "2nne 6 5 4 3 2 1\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    finalize_config(&ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->num_elements, "Number of elements");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->num_nn_levels, "Number of nn levels");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1, se->nn_energy[0], "1st shell nne energy for atom type 0-0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(2, se->nn_energy[1], "1st shell nne energy for atom type 0-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(3, se->nn_energy[2], "1st shell nne energy for atom type 0-2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(4, se->nn_energy[3], "1st shell nne energy for atom type 1-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(5, se->nn_energy[4], "1st shell nne energy for atom type 1-2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(6, se->nn_energy[5], "1st shell nne energy for atom type 2-2");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(6, se->nn_energy[6], "2nd shell nne energy for atom type 0-0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(5, se->nn_energy[7], "2nd shell nne energy for atom type 0-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(4, se->nn_energy[8], "2nd shell nne energy for atom type 0-2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(3, se->nn_energy[9], "2nd shell nne energy for atom type 1-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(2, se->nn_energy[10],
                                     "2nd shell nne energy for atom type 1-2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1, se->nn_energy[11],
                                     "2nd shell nne energy for atom type 2-2");
}

void test_parse_input_cluster_nns_file_success(void) {
    ParseContext ctx = {0};
    FILE *fp = open_file("test/cluster_nns.in");
    TEST_ASSERT_NOT_NULL_MESSAGE(fp, "Failed to open test/cluster_nns.in");

    parse_input_file(fp, &ctx, ss, se, ls);
    finalize_config(&ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_x, "System size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_y, "System size y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, se->system_size_z, "System size z");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2097152, se->max_atoms, "Maximum atom count from system size");

    TEST_ASSERT_EQUAL_INT_MESSAGE(GEOMETRY_CLUSTER, se->geometry, "Geometry type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(32, se->cluster_radius, "Cluster radius");

    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->atom_names_cnt, "Number of atom types");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Ag", se->atom_names[0], "Atom type name at index 0");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Au", se->atom_names[1], "Atom type name at index 1");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, se->substrate_composition[0],
                                     "Composition for atom type 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, se->substrate_composition[1],
                                     "Composition for atom type 1");
    TEST_ASSERT_TRUE_MESSAGE(se->is_soluble[0], "Solubility for atom type 0");
    TEST_ASSERT_FALSE_MESSAGE(se->is_soluble[1], "Solubility for atom type 1");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->dissolution, "Dissolution flag");

    TEST_ASSERT_EQUAL_INT_MESSAGE(FCC, se->lattice_type, "Lattice type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se->num_nn_levels, "Number of nn levels");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[0], "1st shell nne energy at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[1], "1st shell nne energy at index 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.15, se->nn_energy[2], "1st shell nne energy at index 2");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, se->nn_energy[3], "2nd shell nne energy at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, se->nn_energy[4], "2nd shell nne energy at index 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.10, se->nn_energy[5], "2nd shell nne energy at index 2");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.9, se->initial_overpotential, "Initial overpotential");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.03, se->overpotential_ramp_rate,
                                     "Overpotential ramp rate");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.2, se->max_overpotential, "Maximum overpotential");

    TEST_ASSERT_TRUE_MESSAGE(ls->output_iter_csv, "Iter CSV output flag");
    TEST_ASSERT_TRUE_MESSAGE(ls->output_state_csv, "State CSV output flag");
    TEST_ASSERT_TRUE_MESSAGE(ls->output_xyz, "XYZ output flag");
    TEST_ASSERT_EQUAL_INT_MESSAGE(ITERATION_INTERVALS, ls->analysis_type, "Analysis type");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200.0, ls->next_log_checkpoint, "Next log checkpoint");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200.0, ls->log_interval, "Log interval");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->framenum, "Frame number reset");

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293.0, ss->temperature, "Temperature value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(12345U, se->rand_seed, "Random seed value");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FLAVOR_KMC, se->flavor, "Simulation flavor");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, ss->sim_end_type, "Simulation end type");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(2000U, (unsigned int)ss->final_iteration,
                                   "Final iteration value");
}

void test_parse_input_multi_command_unknown_mid_file_fails(void) {
    const char *input = "systemsize 16 16 16\n"
                        "geometry cluster 4\n"
                        "atomtype A B\n"
                        "composition 0.5 0.5\n"
                        "bogus 123\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected unknown command failure in multi-command file");
    });
}

void test_parse_input_multi_command_finalize_composition_sum_fails(void) {
    const char *input = "systemsize 16 16 16\n"
                        "geometry cluster 4\n"
                        "atomtype A B\n"
                        "composition 0.8 0.8\n"
                        "dissolution true false\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_config(&ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected composition sum failure at finalize for multi-command file");
    });
}

void test_parse_input_missing_atomtype_fails(void) {
    const char *input = "composition 1.0\n"
                        "nnlevels 1\n"
                        "1nne 0.1\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_config(&ctx, ss, se, ls);
        fclose(fp);
        TEST_FAIL_MESSAGE("Expected failure");
    });
}

void test_parse_input_missing_nne_shell_fails(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 2\n"
                        "1nne 0.1 0.2 0.3\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_config(&ctx, ss, se, ls);
        fclose(fp);
        TEST_FAIL_MESSAGE("Expected failure");
    });
}

void test_parse_input_unknown_command_fails(void) {
    const char *input = "atomtype A B\n"
                        "foobar 123\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected an unknown command failure");
    });
}

void test_nnlevels_wrong_argcount_fails(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 1 2\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected nnlevels argcount failure");
    });
}

void test_composition_non_numeric_fails(void) {
    const char *input = "atomtype A B\n"
                        "composition 0.5 foo\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected non-numeric composition failure");
    });
}

void test_dissolution_invalid_boolean_fails(void) {
    const char *input = "atomtype A B\n"
                        "dissolution true maybe\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected invalid boolean failure");
    });
}

void test_nne_level_exceeds_nnlevels_fails(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "2nne 0.1 0.2 0.3\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_config(&ctx, ss, se, ls);
        fclose(fp);
        TEST_FAIL_MESSAGE("Expected nne level overflow failure");
    });
}

void test_nne_duplicate_definition_fails(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n"
                        "1nne 0.3 0.2 0.1\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_config(&ctx, ss, se, ls);
        fclose(fp);
        TEST_FAIL_MESSAGE("Expected duplicate nne failure");
    });
}

void test_nne_too_few_values_fails(void) {
    const char *input = "atomtype A B C\n"
                        "nnlevels 1\n"
                        "1nne 1 2 3\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_config(&ctx, ss, se, ls);
        fclose(fp);
        TEST_FAIL_MESSAGE("Expected nne value count failure");
    });
}

void test_composition_count_mismatch_fails(void) {
    const char *input = "atomtype A B C\n"
                        "composition 0.5 0.5\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_config(&ctx, ss, se, ls);
        fclose(fp);
        TEST_FAIL_MESSAGE("Expected composition count mismatch");
    });
}

void test_systemsize_valid_values_success(void) {
    const char *input = "systemsize 10 11 12\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(10, se->system_size_x, "System size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(11, se->system_size_y, "System size y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12, se->system_size_z, "System size z");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1320, se->max_atoms, "Maximum atom count from system size");
}

void test_systemsize_non_numeric_fails(void) {
    const char *input = "systemsize 10 a 12\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected systemsize parse failure");
    });
}

void test_temp_valid_value_success(void) {
    const char *input = "temp 293.5\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(293.5, ss->temperature, "Temperature value");
}

void test_temp_non_numeric_fails(void) {
    const char *input = "temp abc\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected temp parse failure");
    });
}

void test_seed_default_keyword_success(void) {
    const char *input = "seed default\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(DEFAULT_SEED, se->rand_seed, "Default random seed");
}

void test_seed_non_numeric_fails(void) {
    const char *input = "seed notanumber\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected seed parse failure");
    });
}

void test_potential_sweep_values_success(void) {
    const char *input = "potential 0.9 0.03 1.2\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.9, se->initial_overpotential, "Initial overpotential");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.03, se->overpotential_ramp_rate,
                                     "Overpotential ramp rate");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.2, se->max_overpotential, "Maximum overpotential");
}

void test_potential_bad_argcount_fails(void) {
    const char *input = "potential 0.9 0.03\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected potential argcount failure");
    });
}

void test_datalog_interval_values_success(void) {
    const char *input = "datalog iteration interval 200\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(ITERATION_INTERVALS, ls->analysis_type, "Analysis type");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200.0, ls->next_log_checkpoint, "Next log checkpoint");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(200.0, ls->log_interval, "Log interval");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls->framenum, "Frame number reset");
}

void test_datalog_unknown_type_fails(void) {
    const char *input = "datalog weird interval 1\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected datalog type failure");
    });
}

void test_struct_fcc_success(void) {
    const char *input = "struct FCC\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(FCC, se->lattice_type, "Lattice type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(MAXIMUM_NUMBER_OF_NEIGHBORS, se->num_transition_vectors,
                                  "Transition vector count");
}

void test_struct_invalid_value_fails(void) {
    const char *input = "struct HCP\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected struct value failure");
    });
}

void test_output_valid_path_success(void) {
    const char *input = "output output/test_file.out\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_STRING_MESSAGE("output/test_file.out", outFile, "Output filename");
}

void test_output_missing_arg_fails(void) {
    const char *input = "output\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected output arg failure");
    });
}

void test_geometry_cluster_success(void) {
    const char *input = "geometry cluster 32\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(GEOMETRY_CLUSTER, se->geometry, "Geometry type");
    TEST_ASSERT_EQUAL_INT_MESSAGE(32, se->cluster_radius, "Cluster radius");
}

void test_geometry_invalid_type_fails(void) {
    const char *input = "geometry cube 4\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected geometry type failure");
    });
}

void test_run_iteration_mode_success(void) {
    const char *input = "run iteration 2000\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, ss->sim_end_type, "Simulation end type");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(2000U, (unsigned int)ss->final_iteration,
                                   "Final iteration value");
}

void test_run_unknown_mode_fails(void) {
    const char *input = "run steps 100\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected run mode failure");
    });
}

void test_flavor_kmc_success(void) {
    const char *input = "flavor KMC\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(FLAVOR_KMC, se->flavor, "Simulation flavor");
}

void test_flavor_invalid_value_fails(void) {
    const char *input = "flavor BAD\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected flavor value failure");
    });
}

void test_logtype_all_formats_success(void) {
    const char *input = "logtype iter csv xyz\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_TRUE_MESSAGE(ls->output_iter_csv, "Iter CSV output flag");
    TEST_ASSERT_TRUE_MESSAGE(ls->output_state_csv, "State CSV output flag");
    TEST_ASSERT_TRUE_MESSAGE(ls->output_xyz, "XYZ output flag");
}

void test_logtype_missing_arg_fails(void) {
    const char *input = "logtype\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected logtype arg failure");
    });
}

void test_systemsize_wrong_argcount_fails(void) {
    const char *input = "systemsize 10 11\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected systemsize argcount failure");
    });
}

void test_temp_extra_arg_fails(void) {
    const char *input = "temp 300 301\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected temp argcount failure");
    });
}

void test_seed_missing_arg_fails(void) {
    const char *input = "seed\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected seed missing arg failure");
    });
}

void test_potential_non_numeric_sweep_fails(void) {
    const char *input = "potential 0.9 foo 1.2\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected potential numeric parse failure");
    });
}

void test_datalog_invalid_mode_keyword_fails(void) {
    const char *input = "datalog linear cadence 1\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected datalog mode keyword failure");
    });
}

void test_datalog_interval_non_numeric_fails(void) {
    const char *input = "datalog linear interval nope\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected datalog interval parse failure");
    });
}

void test_geometry_cluster_non_numeric_fails(void) {
    const char *input = "geometry cluster radius\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected geometry cluster parse failure");
    });
}

void test_geometry_file_missing_name_fails(void) {
    const char *input = "geometry file\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected geometry file missing name failure");
    });
}

void test_run_iteration_non_numeric_fails(void) {
    const char *input = "run iteration ten\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected run iteration parse failure");
    });
}

void test_run_bad_argcount_fails(void) {
    const char *input = "run time\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected run argcount failure");
    });
}

void test_nne_non_numeric_value_fails(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne 0.1 x 0.3\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected nne numeric parse failure");
    });
}

void test_composition_sum_not_one_fails(void) {
    const char *input = "atomtype A B\n"
                        "composition 0.6 0.5\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_atom_dependent(&ctx, se);
        TEST_FAIL_MESSAGE("Expected composition sum finalize failure");
    });
}

void test_dissolution_count_mismatch_additional_fails(void) {
    const char *input = "atomtype A B C\n"
                        "dissolution true false\n";

    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_atom_dependent(&ctx, se);
        TEST_FAIL_MESSAGE("Expected dissolution count mismatch finalize failure");
    });
}

void test_atomtype_three_types_success(void) {
    const char *input = "atomtype Ag Au Cu\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->atom_names_cnt, "Number of atom names");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->num_elements, "Number of elements");
    TEST_ASSERT_EQUAL_INT_MESSAGE(6, se->num_bond_types, "Number of bond types");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Ag", se->atom_names[0], "Atom type name at index 0");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Au", se->atom_names[1], "Atom type name at index 1");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Cu", se->atom_names[2], "Atom type name at index 2");
}

void test_atomtype_missing_args_fails(void) {
    const char *input = "atomtype\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected atomtype missing-args failure");
    });
}

void test_composition_two_types_finalize_success(void) {
    const char *input = "atomtype A B\n"
                        "composition 0.25 0.75\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    finalize_atom_dependent(&ctx, se);
    fclose(fp);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, se->substrate_composition[0], "Composition at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, se->substrate_composition[1], "Composition at index 1");
}

void test_composition_missing_values_fails(void) {
    const char *input = "composition\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected composition missing-values failure");
    });
}

void test_dissolution_all_true_sets_flag_success(void) {
    const char *input = "atomtype A B C\n"
                        "dissolution true true true\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    finalize_atom_dependent(&ctx, se);
    fclose(fp);

    TEST_ASSERT_TRUE_MESSAGE(se->is_soluble[0], "Solubility at index 0");
    TEST_ASSERT_TRUE_MESSAGE(se->is_soluble[1], "Solubility at index 1");
    TEST_ASSERT_TRUE_MESSAGE(se->is_soluble[2], "Solubility at index 2");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->dissolution, "Dissolution flag");
}

void test_dissolution_missing_values_fails(void) {
    const char *input = "dissolution\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected dissolution missing-values failure");
    });
}

void test_nnlevels_single_value_success(void) {
    const char *input = "nnlevels 3\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->num_nn_levels, "NN levels value");
}

void test_nnlevels_non_numeric_fails(void) {
    const char *input = "nnlevels abc\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected nnlevels non-numeric failure");
    });
}

void test_nne_missing_values_fails(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        TEST_FAIL_MESSAGE("Expected nne missing-values failure");
    });
}

void test_finalize_nne_missing_nnlevels_fails(void) {
    const char *input = "atomtype A B\n"
                        "1nne 0.1 0.2 0.3\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_nne(&ctx, se);
        TEST_FAIL_MESSAGE("Expected finalize_nne without nnlevels failure");
    });
}

void test_finalize_atom_dependent_direct_success(void) {
    const char *input = "atomtype A B\n"
                        "composition 0.4 0.6\n"
                        "dissolution false true\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    finalize_atom_dependent(&ctx, se);
    fclose(fp);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.4, se->substrate_composition[0],
                                     "Finalized composition at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.6, se->substrate_composition[1],
                                     "Finalized composition at index 1");
    TEST_ASSERT_FALSE_MESSAGE(se->is_soluble[0], "Finalized solubility at index 0");
    TEST_ASSERT_TRUE_MESSAGE(se->is_soluble[1], "Finalized solubility at index 1");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, se->dissolution, "Finalized dissolution flag");
}

void test_finalize_atom_dependent_missing_atomtype_fails(void) {
    const char *input = "composition 1.0\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_atom_dependent(&ctx, se);
        TEST_FAIL_MESSAGE("Expected finalize_atom_dependent missing atomtype failure");
    });
}

void test_finalize_nne_direct_single_shell_success(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne 0.11 0.22 0.33\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    parse_input_file(fp, &ctx, ss, se, ls);
    finalize_nne(&ctx, se);
    fclose(fp);

    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se->num_nn_types, "Total NN energy types");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.11, se->nn_energy[0], "NN energy at index 0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.22, se->nn_energy[1], "NN energy at index 1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.33, se->nn_energy[2], "NN energy at index 2");
}

void test_finalize_nne_direct_duplicate_shell_fails(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "1nne 0.1 0.2 0.3\n"
                        "1nne 0.3 0.2 0.1\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_nne(&ctx, se);
        TEST_FAIL_MESSAGE("Expected finalize_nne duplicate shell failure");
    });
}

void test_finalize_nne_direct_level_exceeds_nnlevels_fails(void) {
    const char *input = "atomtype A B\n"
                        "nnlevels 1\n"
                        "2nne 0.1 0.2 0.3\n";
    ParseContext ctx = {0};
    FILE *fp = open_mem(input);

    EXPECT_EXIT(EXIT_FAILURE, {
        parse_input_file(fp, &ctx, ss, se, ls);
        finalize_nne(&ctx, se);
        TEST_FAIL_MESSAGE("Expected finalize_nne level overflow failure");
    });
}


/* ================= Runner ================= */

int main(void) {
    UNITY_BEGIN();

    // End-to-end parser baselines
    RUN_TEST(test_parse_input_two_atomtypes_one_shell_success);
    RUN_TEST(test_parse_input_three_atomtypes_two_shells_success);

    RUN_TEST(test_parse_input_cluster_nns_file_success);
    RUN_TEST(test_parse_input_multi_command_unknown_mid_file_fails);

    RUN_TEST(test_parse_input_multi_command_finalize_composition_sum_fails);
    RUN_TEST(test_parse_input_missing_atomtype_fails);
    RUN_TEST(test_parse_input_missing_nne_shell_fails);
    RUN_TEST(test_parse_input_unknown_command_fails);

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

    // datalog
    RUN_TEST(test_datalog_interval_values_success);
    RUN_TEST(test_datalog_unknown_type_fails);
    RUN_TEST(test_datalog_invalid_mode_keyword_fails);
    RUN_TEST(test_datalog_interval_non_numeric_fails);

    // struct
    RUN_TEST(test_struct_fcc_success);
    RUN_TEST(test_struct_invalid_value_fails);

    // output
    RUN_TEST(test_output_valid_path_success);
    RUN_TEST(test_output_missing_arg_fails);

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

    // logtype
    RUN_TEST(test_logtype_all_formats_success);
    RUN_TEST(test_logtype_missing_arg_fails);

    // Direct finalizer tests
    RUN_TEST(test_finalize_nne_missing_nnlevels_fails);
    RUN_TEST(test_finalize_atom_dependent_direct_success);
    RUN_TEST(test_finalize_atom_dependent_missing_atomtype_fails);
    RUN_TEST(test_finalize_nne_direct_single_shell_success);
    RUN_TEST(test_finalize_nne_direct_duplicate_shell_fails);
    RUN_TEST(test_finalize_nne_direct_level_exceeds_nnlevels_fails);

    int num_failures = UNITY_END();
    
    // return 0 else makefile throws error
    return 0;
}
