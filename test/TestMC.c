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
    // fclose needs to be here in case a test fails
    // set to null to prevent double free from fclose + clean_and_error
    ls->sim_log = NULL;
    fclose(temp_log);
    clean_and_error(0);
}

void test_simulation(void)
{
    ls->sim_log = temp_log;
    char filename[] = "test/mc.in";
    struct SimulationConfig inputs = {0};
    simulation_parameters_from_file(filename, &inputs, ls);
    open_log_files(ls, se->flavor);

    initialize_simulation_from_input(&inputs, ss, se, ls);

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
