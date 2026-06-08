// #include "TTempalte.h"
#include "Mesosim.h"
#include "State.h"
#include "unity.h"
#include <stdlib.h>

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

void test_mesosim(void)
{
    char *argv[] = {"mesosim", "-v", "1000", "test/cluster_nns.in"};
    int argc = sizeof(argv) / sizeof(argv[0]);
    main_mesosim(argc, argv);
    TEST_PASS();
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_mesosim);
    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
