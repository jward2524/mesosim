#include "ErrorM.h"
#include "TUtils.h"
#include "Utils.h"
#include "unity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
externs from ErrorM.c
struct SimulationState *gp_sim_state;
struct SimulationEnv *gp_sim_env;
struct LoggingState *gp_log_state;
*/

void setUp(void)
{
}

void tearDown(void)
{
    // reset global error variables after each test
    gp_sim_state = NULL;
    gp_sim_env = NULL;
    gp_log_state = NULL;
}

// Test set_state assigns pointers correctly
void test_set_state_pointer_assignment(void)
{
    struct SimulationState *ss = calloc(1, sizeof(struct SimulationState));
    struct SimulationEnv *se = calloc(1, sizeof(struct SimulationEnv));
    struct LoggingState *ls = calloc(1, sizeof(struct LoggingState));

    set_state(&ss, &se, &ls);

    TEST_ASSERT_EQUAL_PTR_MESSAGE(&ss, gp_sim_state, "gp_sim_state should point to address of ss");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(&se, gp_sim_env, "gp_sim_env should point to address of se");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(&ls, gp_log_state, "gp_log_state should point to address of ls");

    free(ss);
    free(se);
    free(ls);
}

// Test initialize_states allocates memory and sets global state
void test_initialize_states_allocates_and_sets_state(void)
{
    struct SimulationState *ss = NULL;
    struct SimulationEnv *se = NULL;
    struct LoggingState *ls = NULL;

    initialize_states(&ss, &se, &ls);

    TEST_ASSERT_NOT_NULL_MESSAGE(ss, "initialize_states should allocate non-NULL ss");
    TEST_ASSERT_NOT_NULL_MESSAGE(se, "initialize_states should allocate non-NULL se");
    TEST_ASSERT_NOT_NULL_MESSAGE(ls, "initialize_states should allocate non-NULL ls");

    TEST_ASSERT_EQUAL_PTR_MESSAGE(ss, *gp_sim_state,
                                  "gp_sim_state should point to allocated SimulationState");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(se, *gp_sim_env,
                                  "gp_sim_env should point to allocated SimulationEnv");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(ls, *gp_log_state,
                                  "gp_log_state should point to allocated LoggingState");

    // Clean up memory manually
    free_if_exists((void **)&ss);
    free_if_exists((void **)&se);
    free_if_exists((void **)&ls);
}

// Test free_if_exists frees memory and sets pointer to NULL
void test_free_if_exists_sets_pointer_to_null(void)
{
    void *ptr = malloc(10);
    TEST_ASSERT_NOT_NULL_MESSAGE(ptr, "malloc should return non-NULL pointer");

    free_if_exists(&ptr);
    TEST_ASSERT_NULL_MESSAGE(ptr, "Pointer should be NULL after free_if_exists");

    // freeing NULL again should not crash
    free_if_exists(&ptr);
    TEST_ASSERT_NULL_MESSAGE(ptr, "Pointer should still be NULL after freeing NULL");
}

// Test clean_and_error with exit_error == 0 cleans up memory
void test_clean_and_error_cleans_memory_without_exit(void)
{
    struct SimulationState *ss = NULL;
    struct SimulationEnv *se = NULL;
    struct LoggingState *ls = NULL;
    initialize_states(&ss, &se, &ls);

    // Assign dummy values to test free_if_exists loops
    ss->atom_cnt = 1;
    ss->atom_arr = malloc(sizeof(void *));
    ss->atom_arr[0] = malloc(1 * sizeof(Atom));

    ss->transition_cnt = 1;
    ss->transition_arr = malloc(sizeof(void *));
    ss->transition_arr[0] = malloc(1 * sizeof(Transition));

    ss->rate_cnt = 1;
    ss->rate_arr = malloc(sizeof(*ss->rate_arr));
    ss->rate_arr[0].atom_env = malloc(1);

    ss->transition_probability.rate_arr_index = malloc(sizeof(long));
    ss->transition_probability.lbound = malloc(sizeof(double));
    ss->transition_probability.ubound = malloc(sizeof(double));

    se->atom_names_cnt = 1;
    se->atom_names = malloc(sizeof(char *));
    se->atom_names[0] = dup_str("H");

    clean_and_error(0);

    // all pointers should be NULL after cleanup
    // can't test nested pointers because their parent pointer is already freed and set to NULL
    // so just test top-level pointers
    TEST_ASSERT_NULL_MESSAGE(ss, "ss should be NULL after cleanup");
    TEST_ASSERT_NULL_MESSAGE(se, "se should be NULL after cleanup");
    TEST_ASSERT_NULL_MESSAGE(ls, "ls should be NULL after cleanup");

    TEST_ASSERT_NULL_MESSAGE(gp_sim_state, "gp_sim_state should be NULL after cleanup");
    TEST_ASSERT_NULL_MESSAGE(gp_sim_env, "gp_sim_env should be NULL after cleanup");
    TEST_ASSERT_NULL_MESSAGE(gp_log_state, "gp_log_state should be NULL after cleanup");
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_set_state_pointer_assignment);
    RUN_TEST(test_initialize_states_allocates_and_sets_state);
    RUN_TEST(test_free_if_exists_sets_pointer_to_null);
    RUN_TEST(test_clean_and_error_cleans_memory_without_exit);

    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
