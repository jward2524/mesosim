#include "unity.h"
#include "Input.h"
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

void test_parse_two_atomtypes_one_shell(void) {
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

    TEST_ASSERT_EQUAL_INT(1, se->is_soluble[0]);
    TEST_ASSERT_EQUAL_INT(0, se->is_soluble[1]);

    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.1, se->nn_energy[0], "1nne energy for atom type 0-0");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.2, se->nn_energy[1], "1nne energy for atom type 0-1");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.3, se->nn_energy[2], "1nne energy for atom type 1-1");
}

void test_parse_three_atomtypes_two_shells(void) {
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

void test_missing_atomtype_fails(void) {
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

void test_missing_nne_shell_fails(void) {
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

void test_unknown_command_fails(void) {
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

void test_duplicate_nne_definition_fails(void) {
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


/* ================= Runner ================= */

int main(void) {
    UNITY_BEGIN();

    RUN_TEST(test_parse_two_atomtypes_one_shell);
    RUN_TEST(test_parse_three_atomtypes_two_shells);
    RUN_TEST(test_missing_atomtype_fails);
    RUN_TEST(test_missing_nne_shell_fails);
    RUN_TEST(test_unknown_command_fails);
    RUN_TEST(test_nnlevels_wrong_argcount_fails);
    RUN_TEST(test_composition_non_numeric_fails);
    RUN_TEST(test_dissolution_invalid_boolean_fails);
    RUN_TEST(test_nne_level_exceeds_nnlevels_fails);
    RUN_TEST(test_duplicate_nne_definition_fails);
    RUN_TEST(test_nne_too_few_values_fails);
    RUN_TEST(test_composition_count_mismatch_fails);

    return UNITY_END();
}
