// #include "TTempalte.h"
#include "State.h"
#include "unity.h"

static struct SimulationState *ss;
static struct SimulationEnv *se;
static struct LoggingState *ls;

void setUp(void)
{
    ss = calloc(1, sizeof(struct SimulationState));
    se = calloc(1, sizeof(struct SimulationEnv));
    ls = calloc(1, sizeof(struct LoggingState));
}

void tearDown(void)
{
}

void test_test(void)
{
    TEST_PASS();
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_test);
    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
