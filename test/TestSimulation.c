#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "Simulation.h"
#include "State.h"
#include "TUtils.h"
#include "unity.h"
#include <errno.h>
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
    char filename[] = "test/cluster_nns.in";
    struct SimulationConfig inputs = {0};
    simulation_parameters_from_file(filename, &inputs, ls);
    open_log_files(ls, se->flavor);

    initialize_simulation(&inputs, ss, se, ls);

    perform_kmc(ss, se, ls);
    TEST_PASS();
}

// TODO: add test for number of files produced and number of lines in csv files and values in csv
// files

// TODO: test individual function in Simulation.c
// add_to/remove_from_transition_array cases: removed transition is [only one in list, last in list,
// first in list, middle of list, first in rate list, last in rate list, middle of rate list, only
// one in rate list]

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_simulation);
    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
