// #include "TTempalte.h"
#include "Mesosim.h"
#include "State.h"
#include "TUtils.h"
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

    FILE *f = fopen("test/output/cluster.csv", "rb");
    rewind(f);
    char line[256];
    char *ptr = fgets(line, sizeof(line), f);
    TEST_ASSERT_NOT_NULL_MESSAGE(ptr, "fgets should read a line");
    ptr = fgets(line, sizeof(line), f);
    TEST_ASSERT_NOT_NULL_MESSAGE(ptr, "fgets should read a line");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("0,0,0e+00,226704.000001,293.000000\n", line, "CSV state log");
    fclose(f);

    FILE *checkpoint_file = fopen("test/output/cluster.bin", "rb");
    TEST_ASSERT_NULL_MESSAGE(checkpoint_file, "Checkpoint file should have been deleted");

    const char *output_filenames[] = {
        "test/output/cluster_0_i0.xyz",    "test/output/cluster_1_i500.xyz",
        "test/output/cluster_2_i1000.xyz", "test/output/cluster_3_i1500.xyz",
        "test/output/cluster_4_i2000.xyz", "test/output/cluster_5_i2000.xyz",
        "test/output/cluster.csv",
    };

    assert_many_files_exist_and_remove(output_filenames,
                                       sizeof(output_filenames) / sizeof(output_filenames[0]));
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_mesosim);

    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
