#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "MC.h"
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

void test_mc(void)
{
    ls->sim_log = temp_log;
    char filename[] = "test/mc.in";
    struct SimulationConfig inputs = {0};
    simulation_parameters_from_file(filename, &inputs, ls);
    open_log_files(ls, se->flavor);

    initialize_simulation_from_input(&inputs, ss, se, ls);

    for (int i = 0; i < ss->atom_cnt; ++i) {
        refresh_transitions(i, ss, se);
    }

    perform_metropolis_mc(ss, se, ls);

    rewind(ls->out_formats[0].csv.file);
    char line[256];
    char *ptr = fgets(line, sizeof(line), ls->out_formats[0].csv.file);
    TEST_ASSERT_NOT_NULL_MESSAGE(ptr, "fgets should read a line");
    ptr = fgets(line, sizeof(line), ls->out_formats[0].csv.file);
    TEST_ASSERT_NOT_NULL_MESSAGE(ptr, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("0,1,193941,170451.449999,293.000000\n", line,
                                     "CSV state log after MC step");

    FILE *checkpoint_file = fopen(ls->checkpoint.filename, "rb");
    TEST_ASSERT_NOT_NULL_MESSAGE(checkpoint_file, "Checkpoint file should be created");
    fclose(checkpoint_file);

    fclose(ls->out_formats[0].csv.file);
    ls->out_formats[0].csv.file = NULL;
    int ret = remove(ls->out_formats[0].csv.filename);
    if (ret != 0) {
        fprintf(stderr, "Error deleting file %s: %s\n", ls->out_formats[0].csv.filename,
                strerror(errno));
    }
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_mc);
    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
