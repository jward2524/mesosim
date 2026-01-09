#include "MC.h"
#include "FileIO.h"
#include "State.h"
#include "unity.h"
#include "ErrorM.h"
#include "Initialization.h"
#include <string.h>
#include <stdlib.h>

struct SimulationState *ss;
struct SimulationEnv *se;
struct LoggingState *ls;
FILE *temp_log;
const char temp_name[] = "temp.log";

static void fopen_error(const char *filename, FILE *file)
{
    if (file == NULL) {
        fprintf(stderr, "Couldn't open file %s: %s\n", filename, strerror(errno));
        TEST_ASSERT_NOT_NULL_MESSAGE(file, "File not opened - check result file");
    }
}

void setUp(void)
{
    ss = calloc(1, sizeof(struct SimulationState));
    se = calloc(1, sizeof(struct SimulationEnv));
    ls = calloc(1, sizeof(struct LoggingState));
    
    set_state(ss, se, ls);
    temp_log = fopen(temp_name, "w");
    fopen_error(temp_name, temp_log);
}

void tearDown(void) {
    free(ss);
    free(se);
    free(ls);

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
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_simulation);
    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
