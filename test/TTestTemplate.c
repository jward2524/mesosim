// #include "TTempalte.h"
#include "State.h"
#include "unity.h"

struct SimulationState *ss;
struct SimulationEnv *se;
struct LoggingState *ls;

void setUp(void)
{
    ss = calloc(1, sizeof(struct SimulationState));
    se = calloc(1, sizeof(struct SimulationEnv));
    ls = calloc(1, sizeof(struct LoggingState));
}

void tearDown(void) {}

void test_parse_comment(void)
{
    char line[] = "systemsize 128 128 128";

    TEST_ASSERT_EQUAL_INT_MESSAGE(128, 128, "system size x");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, 128, "system size y");
    TEST_ASSERT_EQUAL_INT_MESSAGE(128, 128, "system size z");
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_parse_comment);
    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
