#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "MC.h"
#include "State.h"
#include "TUtils.h"
#include "unity.h"
#include <stdlib.h>
#include <string.h>

struct SimulationState *ss;
struct SimulationEnv *se;
struct LoggingState *ls;
FILE *temp_log;

void setUp(void)
{
    initialize_states(&ss, &se, &ls);
    init_temp(&temp_log);
}

void tearDown(void)
{
    clean_and_error(0);

    // fclose needs to be here in case a test fails
    fclose(temp_log);
}

void test_simulation(void)
{
    char filename[] = "test/mc.in";
    time_t starttime = 0;
    time(&starttime);
    simulation_parameters_from_file(filename, ss, se, ls, temp_log, starttime);

    initialize_simulation(ss, se, ls);

    perform_metropolis_mc(ss, se, ls);
    TEST_PASS();
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_simulation);
    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
