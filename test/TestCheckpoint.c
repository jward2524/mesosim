#include "Checkpoint.h"
#include "Initialization.h"
#include "TUtils.h"
#include "unity.h"
#include <stdlib.h>
#include <string.h>

FILE *temp_checkpoint_file;
const char temp_checkpoint_path[] = "build/test/temp_checkpoint.bin";

void setUp(void)
{
}

void tearDown(void)
{
    int pexists = temp_checkpoint_file != NULL;
    close_if_exists(&temp_checkpoint_file);
    int rc = remove(temp_checkpoint_path);

    // if file pointer didn't exist, file still could have existed bc created elsewhere
    // if file pointer existed and remove failed, print error
    // if file pointer didn't exist, file still may or may not have existed
    // if file didn't exist, remove will fail, but is ok
    if (pexists && rc) {
        perror("Remove of test checkpoint file failed");
    }
}

void test_checkpoint_checksum32_is_deterministic(void)
{
    const unsigned char payload[] = {0x01u, 0x02u, 0x03u, 0x04u};
    const unsigned char payload2[] = {0x01u, 0x02u, 0x03u, 0x05u};
    uint32_t checksum1 = checkpoint_checksum32(payload, sizeof(payload));
    uint32_t checksum2 = checkpoint_checksum32(payload, sizeof(payload));
    uint32_t checksum3 = checkpoint_checksum32(payload2, sizeof(payload2));

    TEST_ASSERT_EQUAL_UINT32_MESSAGE(checksum1, checksum2,
                                     "checksum should be stable for the same bytes");
    TEST_ASSERT_NOT_EQUAL_UINT32_MESSAGE(checksum1, checksum3,
                                         "checksum should change when the payload changes");
}

void test_checkpoint_header_magic_helpers(void)
{
    CheckpointHeader header = {0};

    checkpoint_header_set_magic(&header);

    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_has_valid_magic(&header),
                             "magic should match after being set");
    header.magic[0] = 'X';
    TEST_ASSERT_FALSE_MESSAGE(checkpoint_header_has_valid_magic(&header),
                              "magic should fail after corruption");
}

void test_checkpoint_header_init_and_finalize(void)
{
    CheckpointHeader header;

    checkpoint_header_init(&header);
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(CHECKPOINT_FORMAT_VERSION, header.format_version,
                                     "header version should be initialized");
    TEST_ASSERT_EQUAL_UINT32_MESSAGE((uint32_t)sizeof(CheckpointHeader), header.header_bytes,
                                     "header size should be stored");
    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_has_valid_structure(&header),
                             "header shape should be valid after init");

    checkpoint_header_finalize(&header, 128u, 0x12345678u);
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(128u, header.payload_bytes, "payload size should be written");
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(0x12345678u, header.checksum, "checksum should be written");
}

void test_checkpoint_header_write_and_read_round_trip(void)
{
    CheckpointHeader written;
    CheckpointHeader read_back;

    temp_checkpoint_file = fopen(temp_checkpoint_path, "w+");
    fopen_error(temp_checkpoint_path, temp_checkpoint_file);

    checkpoint_header_init(&written);
    checkpoint_header_finalize(&written, 256u, 0xDEADBEEFu);

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK,
                                  checkpoint_header_write(temp_checkpoint_file, &written),
                                  "header should write successfully");
    rewind(temp_checkpoint_file);

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK,
                                  checkpoint_header_read(temp_checkpoint_file, &read_back),
                                  "header should read successfully");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&written, &read_back, sizeof(CheckpointHeader),
                                     "header should round-trip byte-for-byte");
}

void test_checkpoint_save_and_load_header_only(void)
{
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    save_status = checkpoint_save(temp_checkpoint_path, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "header-only checkpoint should save successfully");

    load_status = checkpoint_load(temp_checkpoint_path, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "header-only checkpoint should load successfully");
}

void test_checkpoint_save_and_load_state_scalars(void)
{
    struct SimulationState ss = {0};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    ss.iter = 42ul;
    ss.mmc_steps = 99ul;
    ss.final_iteration = 1234ul;
    ss.run_stime = 12.5;
    ss.simulation_should_kill_itself = true;
    ss.elapsed_stime = 3.75;
    ss.sim_end_type = SIM_END_BY_ITERATIONS;
    ss.frequency_sum = 0.125;
    ss.total_internal_energy = -4.5;
    ss.temperature = 295.0;
    ss.overpotential = 0.8;
    ss.total_atoms_dissolved = 7;

    save_status = checkpoint_save(temp_checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "state scalar checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));

    load_status = checkpoint_load(temp_checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "state scalar checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(42u, ss.iter, "iteration counter should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(99u, ss.mmc_steps, "mmc steps should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(1234u, ss.final_iteration, "final iteration should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(12.5, ss.run_stime, "runtime limit should restore");
    TEST_ASSERT_TRUE_MESSAGE(ss.simulation_should_kill_itself, "kill-itself flag should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(3.75, ss.elapsed_stime, "elapsed time should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, ss.sim_end_type,
                                  "simulation end type should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.125, ss.frequency_sum, "frequency sum should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(-4.5, ss.total_internal_energy,
                                     "internal energy should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(295.0, ss.temperature, "temperature should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.8, ss.overpotential, "overpotential should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(7, ss.total_atoms_dissolved,
                                  "dissolved atom count should restore");
}

void test_checkpoint_save_and_load_env_scalars(void)
{
    struct SimulationEnv se = {0};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    // populate a minimal env with key config fields for checkpointing.
    se.flavor = FLAVOR_KMC;
    se.lattice_type = FCC;
    se.rand_seed = 12345u;
    se.system_size_x = 20;
    se.system_size_y = 20;
    se.system_size_z = 20;
    se.num_elements = 2;
    se.num_nn_levels = 2;
    se.num_bond_types = 3;
    se.num_nn_types = 6;
    se.num_transition_vectors = 12;
    se.overpotential_ramp_rate = 0.01;
    se.max_overpotential = 1.0;

    save_status = checkpoint_save(temp_checkpoint_path, NULL, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "env scalar checkpoint should save successfully");

    memset(&se, 0, sizeof(se));

    load_status = checkpoint_load(temp_checkpoint_path, NULL, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "env scalar checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, se.flavor, "flavor should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(12345u, se.rand_seed, "random seed should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(20, se.system_size_x, "system size x should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(20, se.system_size_y, "system size y should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(20, se.system_size_z, "system size z should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se.num_elements, "num elements should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se.num_nn_levels, "num nn levels should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se.num_bond_types, "num bond types should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(6, se.num_nn_types, "num nieghbor types should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12, se.num_transition_vectors,
                                  "num transition vectors should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.01, se.overpotential_ramp_rate,
                                     "overpotential ramp rate should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.0, se.max_overpotential, "max overpotential should restore");

    remove(temp_checkpoint_path);
}

void test_checkpoint_save_and_load_state_and_env(void)
{
    struct SimulationState ss = {0};
    struct SimulationEnv se = {0};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    // populate both state and env for a combined checkpoint
    ss.iter = 100ul;
    ss.mmc_steps = 99ul;
    ss.final_iteration = 1234ul;
    ss.run_stime = 12.5;
    ss.simulation_should_kill_itself = true;
    ss.elapsed_stime = 3.75;
    ss.sim_end_type = SIM_END_BY_ITERATIONS;
    ss.frequency_sum = 0.125;
    ss.total_internal_energy = -4.5;
    ss.temperature = 300.0;
    ss.overpotential = 0.8;
    ss.total_atoms_dissolved = 7;

    se.flavor = FLAVOR_MC;
    se.lattice_type = FCC;
    se.rand_seed = 12345u;
    se.system_size_x = 40;
    se.system_size_y = 20;
    se.system_size_z = 20;
    se.num_elements = 2;
    se.num_nn_levels = 1;
    se.num_bond_types = 3;
    se.num_transition_vectors = 12;
    se.overpotential_ramp_rate = 0.01;
    se.max_overpotential = 1.0;

    save_status = checkpoint_save(temp_checkpoint_path, &ss, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "combined state+env checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));
    memset(&se, 0, sizeof(se));

    load_status = checkpoint_load(temp_checkpoint_path, &ss, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "combined state+env checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(100u, ss.iter, "iter should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(300.0, ss.temperature,
                                     "temperature should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_MC, se.flavor,
                                   "flavor should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(40, se.system_size_x,
                                  "system size should restore from combined checkpoint");
}

void test_checkpoint_save_and_load_atom_round_trip(void)
{
    struct SimulationState ss = {0};
    struct SimulationEnv se = {0};
    Atom original_atom = {.lattice = {1, 2, 3}, .type = 0};
    Atom *atom_refs[1] = {&original_atom};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    // populate both state and env for a combined checkpoint
    ss.iter = 42ul;
    ss.mmc_steps = 99ul;
    ss.final_iteration = 1234ul;
    ss.run_stime = 12.5;
    ss.simulation_should_kill_itself = true;
    ss.elapsed_stime = 3.75;
    ss.sim_end_type = SIM_END_BY_ITERATIONS;
    ss.frequency_sum = 0.125;
    ss.total_internal_energy = -4.5;
    ss.temperature = 295.0;
    ss.overpotential = 0.8;
    ss.total_atoms_dissolved = 7;

    se.flavor = FLAVOR_KMC;
    se.lattice_type = FCC;
    se.rand_seed = 12345u;
    se.system_size_x = 20;
    se.system_size_y = 20;
    se.system_size_z = 20;
    se.zone_count_u = 20;
    se.zone_count_v = 20;
    se.zone_count_w = 20;
    se.num_elements = 2;
    se.num_nn_levels = 1;
    se.num_bond_types = 3;
    se.num_transition_vectors = 12;
    se.overpotential_ramp_rate = 0.01;
    se.max_overpotential = 1.0;
    se.is_soluble = (bool *)malloc((size_t)se.num_elements * sizeof(bool));
    TEST_ASSERT_NOT_NULL_MESSAGE(se.is_soluble, "is_soluble array should allocate");
    for (int i = 0; i < se.num_elements; ++i) {
        se.is_soluble[i] = false;
    }

    ss.atom_cnt = 1;
    ss.atom_arr = atom_refs;

    save_status = checkpoint_save(temp_checkpoint_path, &ss, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "atom checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));

    load_status = checkpoint_load(temp_checkpoint_path, &ss, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "atom checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(1u, (unsigned int)ss.atom_cnt, "atom count should restore");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss.atom_arr, "atom array should be allocated on restore");
    TEST_ASSERT_EQUAL_INT_ARRAY_MESSAGE(original_atom.lattice, ss.atom_arr[0]->lattice, 3,
                                        "atom should round-trip lattice field");
    TEST_ASSERT_EQUAL_INT_MESSAGE(original_atom.type, ss.atom_arr[0]->type,
                                  "atom should round-trip type field");

    free(ss.atom_arr[0]);
    free(ss.atom_arr);
}

void test_checkpoint_load_missing_file_returns_error(void)
{
    CheckpointStatus load_status;

    remove(temp_checkpoint_path);

    load_status = checkpoint_load(temp_checkpoint_path, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, load_status,
                                  "loading a missing checkpoint should fail cleanly");
}

void test_checkpoint_load_corrupted_payload_returns_error(void)
{
    struct SimulationState ss = {0};
    FILE *file;
    int corrupted_byte;
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    ss.iter = 11ul;
    ss.temperature = 273.15;
    ss.overpotential = 0.25;

    save_status = checkpoint_save(temp_checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "baseline checkpoint should save successfully");

    file = fopen(temp_checkpoint_path, "rb+");
    TEST_ASSERT_NOT_NULL_MESSAGE(file, "checkpoint file should open for corruption");

    // move file position to start of file
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, fseek(file, (long)sizeof(CheckpointHeader), SEEK_SET),
                                  "file should seek to payload start");

    // rewrite the first byte of the file with a corrupted byte
    corrupted_byte = fgetc(file);
    TEST_ASSERT_TRUE_MESSAGE(corrupted_byte != EOF, "payload byte should be readable");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, fseek(file, (long)sizeof(CheckpointHeader), SEEK_SET),
                                  "file should seek back to payload start for overwrite");
    // flip bits and write byte to start of file
    TEST_ASSERT_NOT_EQUAL_MESSAGE(EOF, fputc(corrupted_byte ^ 0xFF, file),
                                  "payload byte should be overwritten");
    fclose(file);

    load_status = checkpoint_load(temp_checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, load_status,
                                  "loading a corrupted checkpoint should fail cleanly");
}

void test_checkpoint_rebuild_zones_and_rates_from_atoms(void)
{
    struct SimulationState ss_orig = {0};
    struct SimulationState ss_restored = {0};
    struct SimulationEnv se = {0};
    CheckpointStatus save_status, load_status;

    // setup minimal env to support zone queries
    struct UserInputs dummy_config = {
        .system_size_x = 40,
        .system_size_y = 40,
        .system_size_z = 40,
        .lattice_type = FCC,
        .flavor = FLAVOR_KMC,
        .num_nn_levels = 1,
        .num_elements = 1,
        .num_nn_types = 1,
    };
    dummy_config.is_soluble = (bool *)malloc((size_t)dummy_config.num_elements * sizeof(bool));
    TEST_ASSERT_NOT_NULL_MESSAGE(dummy_config.is_soluble, "is_soluble array should allocate");
    for (int i = 0; i < dummy_config.num_elements; ++i) {
        dummy_config.is_soluble[i] = false;
    }
    dummy_config.nn_energy = (double *)malloc(1 * sizeof(double));
    TEST_ASSERT_NOT_NULL_MESSAGE(dummy_config.nn_energy, "nn_energy array should allocate");
    dummy_config.nn_energy[0] = 0.5;

    se.zone_count_u = 10;
    se.zone_count_v = 10;
    se.zone_count_w = 10;
    initialize_env_from_config(&dummy_config, &se);

    // setup original state with a few atoms
    ss_orig.iter = 50ul;
    ss_orig.temperature = 300.0;
    ss_orig.atom_cnt = 2u;

    // allocate and setup atoms
    ss_orig.atom_arr = (Atom **)malloc(2 * sizeof(Atom *));
    TEST_ASSERT_NOT_NULL(ss_orig.atom_arr);

    for (int i = 0; i < 2; ++i) {
        ss_orig.atom_arr[i] = (Atom *)malloc(sizeof(Atom));
        TEST_ASSERT_NOT_NULL(ss_orig.atom_arr[i]);
        memset(ss_orig.atom_arr[i], 0, sizeof(Atom));

        ss_orig.atom_arr[i]->type = 0;
        ss_orig.atom_arr[i]->energy = -0.5;
        ss_orig.atom_arr[i]->lattice[0] = i;
        ss_orig.atom_arr[i]->lattice[1] = 0;
        ss_orig.atom_arr[i]->lattice[2] = 0;
        ss_orig.atom_arr[i]->next_atom = -1;
        ss_orig.atom_arr[i]->previous_atom = -1;

        // mark all neighbor slots as empty
        for (int j = 0; j < MAXIMUM_NUMBER_OF_NEIGHBORS; ++j) {
            ss_orig.atom_arr[i]->neighbor_atom_idxs[j] = -1;
        }
        for (int j = 0; j < MAXIMUM_NUMBER_OF_NEIGHBORS + DISSOLUTION; ++j) {
            ss_orig.atom_arr[i]->transition_indices[j] = -1;
        }
    }

    // save checkpoint with atoms
    save_status = checkpoint_save(temp_checkpoint_path, &ss_orig, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "checkpoint with atoms should save successfully");

    // allocate env for restore
    initialize_zones(&ss_restored.zone_arr, &se);
    TEST_ASSERT_NOT_NULL(ss_restored.zone_arr);

    // load checkpoint
    load_status = checkpoint_load(temp_checkpoint_path, &ss_restored, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "checkpoint with atoms should load successfully");

    // verify restored state
    TEST_ASSERT_EQUAL_INT_MESSAGE(ss_orig.iter, ss_restored.iter, "iter should be restored");
    TEST_ASSERT_EQUAL_INT_MESSAGE(ss_orig.atom_cnt, ss_restored.atom_cnt,
                                  "atom count should be restored");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss_restored.atom_arr, "atom_arr should be allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss_restored.zone_arr, "zone_arr should be allocated");

    // verify atoms have proper initialization for their transition indices (zeroed during
    // restore)
    for (int i = 0; i < ss_restored.atom_cnt; ++i) {
        TEST_ASSERT_NOT_NULL_MESSAGE(ss_restored.atom_arr[i], "each atom should be allocated");
        // transition_indices should be set during rebuild (at least first one should be
        // attempted)
        TEST_ASSERT_TRUE_MESSAGE(ss_restored.atom_arr[i]->lattice[0] >= 0,
                                 "atom lattice coordinates should be valid");
    }

    // cleanup just the loaded atoms (don't free zone_arr from checkpoint reload)
    for (int i = 0; i < ss_orig.atom_cnt; ++i) {
        free(ss_orig.atom_arr[i]);
    }
    free(ss_orig.atom_arr);

    for (int i = 0; i < ss_restored.atom_cnt; ++i) {
        free(ss_restored.atom_arr[i]);
    }
    free(ss_restored.atom_arr);

    // cleanup zones (these were allocated during restore)
    if (ss_restored.zone_arr) {
        for (size_t i = 0; i < se.zone_count_u; ++i) {
            if (ss_restored.zone_arr[i]) {
                for (size_t j = 0; j < se.zone_count_v; ++j) {
                    free(ss_restored.zone_arr[i][j]);
                }
                free(ss_restored.zone_arr[i]);
            }
        }
        free(ss_restored.zone_arr);
    }
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_checkpoint_checksum32_is_deterministic);
    RUN_TEST(test_checkpoint_header_magic_helpers);
    RUN_TEST(test_checkpoint_header_init_and_finalize);
    RUN_TEST(test_checkpoint_header_write_and_read_round_trip);
    RUN_TEST(test_checkpoint_save_and_load_header_only);
    RUN_TEST(test_checkpoint_save_and_load_state_scalars);
    RUN_TEST(test_checkpoint_save_and_load_env_scalars);
    RUN_TEST(test_checkpoint_save_and_load_state_and_env);
    RUN_TEST(test_checkpoint_save_and_load_atom_round_trip);
    RUN_TEST(test_checkpoint_load_missing_file_returns_error);
    RUN_TEST(test_checkpoint_load_corrupted_payload_returns_error);
    RUN_TEST(test_checkpoint_rebuild_zones_and_rates_from_atoms);

    UNITY_END();
    return 0;
}
