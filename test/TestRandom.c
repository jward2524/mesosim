#include "Random.h"
#include "unity.h"

void setUp(void)
{
}

void tearDown(void)
{
}

void test_sran_initializes_state_deterministically(void)
{
    RandomState state1;
    RandomState state2;

    sran(12345ULL, &state1);
    sran(12345ULL, &state2);

    TEST_ASSERT_EQUAL_UINT64_MESSAGE(state1.u, state2.u,
                                     "sran should produce the same state for the same seed");
    TEST_ASSERT_EQUAL_UINT64_MESSAGE(state1.v, state2.v,
                                     "sran should produce the same state for the same seed");
    TEST_ASSERT_EQUAL_UINT64_MESSAGE(state1.w, state2.w,
                                     "sran should produce the same state for the same seed");

    for (int i = 0; i < 10; i++) {
        unsigned long long r1 = ran(&state1);
        unsigned long long r2 = ran(&state2);

        char msg[256];
        snprintf(msg, sizeof(msg), "ran should produce the same value for the same seed on call %d",
                 i);
        TEST_ASSERT_EQUAL_UINT64_MESSAGE(r1, r2, msg);
    }
}

void test_ran_produces_different_values_on_consecutive_calls(void)
{
    RandomState state;
    sran(12345ULL, &state);

    unsigned long long r1 = ran(&state);
    unsigned long long r2 = ran(&state);

    TEST_ASSERT_NOT_EQUAL_UINT64_MESSAGE(
        r1, r2, "ran should produce different values on consecutive calls");
}

void test_ran_produces_different_values_for_different_seeds(void)
{
    RandomState state1;
    RandomState state2;

    sran(12345ULL, &state1);
    sran(54321ULL, &state2);

    unsigned long long r1 = ran(&state1);
    unsigned long long r2 = ran(&state2);

    TEST_ASSERT_NOT_EQUAL_UINT64_MESSAGE(
        r1, r2, "ran should produce different values for different seeds on the first call");
}

void test_dran_produces_floating_point_values_in_range(void)
{
    RandomState state;
    sran(12345ULL, &state);

    for (int i = 0; i < 1000; i++) {
        double r = dran(&state);
        TEST_ASSERT_TRUE_MESSAGE(r >= 0.0 && r < 1.0,
                                 "dran should produce values in the range [0.0, 1.0)");
    }
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_sran_initializes_state_deterministically);
    RUN_TEST(test_ran_produces_different_values_on_consecutive_calls);
    RUN_TEST(test_ran_produces_different_values_for_different_seeds);
    RUN_TEST(test_dran_produces_floating_point_values_in_range);

    UNITY_END();
    return 0;
}
