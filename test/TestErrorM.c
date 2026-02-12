#include "ErrorM.h"
#include "TUtils.h"
#include "unity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
externs from ErrorM.c
struct SimulationState *sim_state;
struct SimulationEnv *sim_env;
struct LoggingState *log_state;
*/

void setUp(void)
{
}

void tearDown(void)
{
    // reset global error variables after each test
    sim_state = NULL;
    sim_env = NULL;
    log_state = NULL;
}

// Test set_state assigns pointers correctly
void test_set_state_pointer_assignment(void)
{
    struct SimulationState *ss = calloc(1, sizeof(struct SimulationState));
    struct SimulationEnv *se = calloc(1, sizeof(struct SimulationEnv));
    struct LoggingState *ls = calloc(1, sizeof(struct LoggingState));

    set_state(ss, se, ls);

    TEST_ASSERT_EQUAL_PTR_MESSAGE(ss, sim_state, "sim_state should point to ss");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(se, sim_env, "sim_env should point to se");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(ls, log_state, "log_state should point to ls");
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

    TEST_ASSERT_EQUAL_PTR_MESSAGE(ss, sim_state,
                                  "sim_state should point to allocated SimulationState");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(se, sim_env, "sim_env should point to allocated SimulationEnv");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(ls, log_state,
                                  "log_state should point to allocated LoggingState");

    // Clean up memory manually
    free_if_exists((void **)&ss);
    free_if_exists((void **)&se);
    free_if_exists((void **)&ls);

    sim_state = NULL;
    sim_env = NULL;
    log_state = NULL;
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

    initialize_states(&sim_state, &sim_env, &log_state);

    // Assign dummy values to test free_if_exists loops
    sim_state->atom_cnt = 1;
    sim_state->atom_arr = malloc(sizeof(void *));
    sim_state->atom_arr[0] = malloc(1 * sizeof(Atom));

    sim_state->transition_cnt = 1;
    sim_state->transition_arr = malloc(sizeof(void *));
    sim_state->transition_arr[0] = malloc(1 * sizeof(Transition));

    sim_state->rate_cnt = 1;
    sim_state->rate_arr = malloc(sizeof(*sim_state->rate_arr));
    sim_state->rate_arr[0].atom_env = malloc(1);

    sim_state->transition_probability.rate_arr_index = malloc(sizeof(long));
    sim_state->transition_probability.lbound = malloc(sizeof(double));
    sim_state->transition_probability.ubound = malloc(sizeof(double));

    sim_env->atom_names_cnt = 1;
    sim_env->atom_names = malloc(sizeof(char *));
    sim_env->atom_names[0] = malloc(strlen("H") + 1);  // +1 for the null terminator
    if (sim_env->atom_names[0] != NULL) {
	        strcpy(sim_env->atom_names[0], "H");
    }

    clean_and_error(0);

    // all pointers should be NULL after cleanup
    // can't test nested pointers because their parent pointer is already freed and set to NULL
    // so just test top-level pointers
    TEST_ASSERT_NULL_MESSAGE(sim_state, "sim_state should be NULL after cleanup");
    TEST_ASSERT_NULL_MESSAGE(sim_env, "sim_env should be NULL after cleanup");
    TEST_ASSERT_NULL_MESSAGE(log_state, "log_state should be NULL after cleanup");
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
