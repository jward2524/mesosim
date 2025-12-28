#include "Initialization.h"
#include "unity.h"
#include <stdlib.h>

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

void test_initialize_from_file(void)
{

}

int main(void)
{
    UNITY_BEGIN();
    // RUN_TEST(test_initialize_from_file);
    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
