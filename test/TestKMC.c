#include "Checkpoint.h"
#include "ErrorM.h"
#include "FileIO.h"
#include "Initialization.h"
#include "KMC.h"
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
    char filename[] = "test/cluster_nns.in";
    struct SimulationConfig inputs = {0};
    simulation_parameters_from_file(filename, &inputs, ls);
    open_log_files(ls, se->flavor);

    initialize_simulation_from_input(&inputs, ss, se, ls);

    for (int i = 0; i < ss->atom_cnt; ++i) {
        refresh_transitions(i, ss, se);
    }
    compute_transition_array(ss, se);

    perform_kmc(ss, se, ls);

    rewind(ls->out_formats[0].csv.file);
    char line[256];
    char *ptr = fgets(line, sizeof(line), ls->out_formats[0].csv.file);
    TEST_ASSERT_NOT_NULL_MESSAGE(ptr, "fgets should read a line");
    ptr = fgets(line, sizeof(line), ls->out_formats[0].csv.file);
    TEST_ASSERT_NOT_NULL_MESSAGE(ptr, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("0,200,1.01766e+00,226468.050001,293.000000\n", line,
                                     "CSV state log after KMC step");

    fclose(ls->out_formats[0].csv.file);
    ls->out_formats[0].csv.file = NULL;

    FILE *checkpoint_file = fopen(ls->checkpoint.filename, "rb");
    TEST_ASSERT_NOT_NULL_MESSAGE(checkpoint_file, "Checkpoint file should be created");
    fclose(checkpoint_file);

    // validate checkpoint file is valid
    CheckpointStatus status = read_checkpoint_file(ls->checkpoint.filename, ss, se, ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "Checkpoint file should be valid");

    // no initial output '_0_i0' file because that is done by Mesosim.c
    const char *output_filenames[] = {
        "test/output/cluster_0_i500.xyz",  "test/output/cluster_1_i1000.xyz",
        "test/output/cluster_2_i1500.xyz", "test/output/cluster_3_i2000.xyz",
        "test/output/cluster.csv",
    };

    assert_many_files_exist_and_remove(output_filenames,
                                       sizeof(output_filenames) / sizeof(output_filenames[0]));
}

// TODO: test individual function in Simulation.c
// add_to/remove_from_transition_array cases: removed transition is [only one in list, last in list,
// first in list, middle of list, first in rate list, last in rate list, middle of rate list, only
// one in rate list]

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_mc);
    UNITY_END();

    clean_temp(&temp_log);

    // return 0 else makefile throws error
    return 0;
}
