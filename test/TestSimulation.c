#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "Simulation.h"
#include "State.h"
#include "unity.h"
#include <stdlib.h>
#include <string.h>
#include <errno.h>

struct SimulationState *ss;
struct SimulationEnv *se;
struct LoggingState *ls;
FILE *temp_log;
const char temp_name[] = "temp.log";

static void fopen_error(const char *filename, FILE *file)
{
    if (file == NULL) {
        printf("ERROR! Couldn't open output file %s\n", filename);
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

void tearDown(void)
{
    clean_and_exit(0);

    // fclose needs to be here in case a test fails
    fclose(temp_log);
}

void test_simulation(void)
{
    char filename[] = "test/cluster_nns.in";
    time_t starttime = 0;
    time(&starttime);
    simulation_parameters_from_file(filename, ss, se, ls, temp_log, starttime);

    get_shifts(se);
    initialize_zones(ss->zone_arr, se);
    set_primitive_basis(se);
    initialize_simulation_box(se);
    initialize_neighbor_offsets(se);
    initialize_simulation_variables(ss, se);
    initialize_initial_structure(ss, se, ls);

    perform_simulation(ss, se, ls);
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_simulation);
    UNITY_END();
    
    if (temp_log) {
        int rc = remove(temp_name);
        if (rc)
            perror("Remove of test log file failed");
    }

    // return 0 else makefile throws error
    return 0;
}
