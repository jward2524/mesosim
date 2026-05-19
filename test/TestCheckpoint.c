#include "Checkpoint.h"
#include "TUtils.h"
#include "unity.h"
#include <stdlib.h>
#include <string.h>

FILE *temp_log;

void setUp(void)
{
    init_temp(&temp_log);
}

void tearDown(void)
{
    close_if_exists(&temp_log);
    clean_temp(&temp_log);
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
    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_has_valid_shape(&header),
                             "header shape should be valid after init");

    checkpoint_header_finalize(&header, 128u, 0x12345678u);
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(128u, header.payload_bytes, "payload size should be written");
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(0x12345678u, header.checksum, "checksum should be written");
}

void test_checkpoint_header_write_and_read_round_trip(void)
{
    CheckpointHeader written;
    CheckpointHeader read_back;

    checkpoint_header_init(&written);
    checkpoint_header_finalize(&written, 256u, 0xDEADBEEFu);

    rewind(temp_log);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, checkpoint_header_write(temp_log, &written),
                                  "header should write successfully");
    rewind(temp_log);

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, checkpoint_header_read(temp_log, &read_back),
                                  "header should read successfully");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&written, &read_back, sizeof(CheckpointHeader),
                                     "header should round-trip byte-for-byte");
}

void test_checkpoint_save_and_load_header_only(void)
{
    const char checkpoint_path[] = "build/test/checkpoint_header_only.bin";
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    save_status = checkpoint_save(checkpoint_path, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "header-only checkpoint should save successfully");

    load_status = checkpoint_load(checkpoint_path, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "header-only checkpoint should load successfully");

    remove(checkpoint_path);
}

void test_checkpoint_save_and_load_state_scalars(void)
{
    const char checkpoint_path[] = "build/test/checkpoint_state_scalars.bin";
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

    save_status = checkpoint_save(checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "state scalar checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));

    load_status = checkpoint_load(checkpoint_path, &ss, NULL);
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

    remove(checkpoint_path);
}

void test_checkpoint_save_and_load_env_scalars(void)
{
    const char checkpoint_path[] = "build/test/checkpoint_env_scalars.bin";
    struct SimulationEnv se = {0};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    // populate a minimal env with key config fields for checkpointing.
    se.flavor = FLAVOR_KMC;
    se.rand_seed = 12345u;
    se.system_size_x = 20;
    se.system_size_y = 20;
    se.system_size_z = 20;
    se.num_elements = 2;
    se.num_nn_levels = 3;
    se.num_bond_types = 6;
    se.num_neighbor_types = 8;
    se.num_transition_vectors = 12;
    se.primitive_basis[0][0] = 1.0;
    se.primitive_basis[0][1] = 0.0;
    se.primitive_basis[0][2] = 0.0;
    se.primitive_basis[1][0] = 0.0;
    se.primitive_basis[1][1] = 1.0;
    se.primitive_basis[1][2] = 0.0;
    se.primitive_basis[2][0] = 0.0;
    se.primitive_basis[2][1] = 0.0;
    se.primitive_basis[2][2] = 1.0;
    se.overpotential_ramp_rate = 0.01;
    se.max_overpotential = 1.0;

    save_status = checkpoint_save(checkpoint_path, NULL, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "env scalar checkpoint should save successfully");

    memset(&se, 0, sizeof(se));

    load_status = checkpoint_load(checkpoint_path, NULL, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "env scalar checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, se.flavor, "flavor should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(12345u, se.rand_seed, "random seed should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(20, se.system_size_x, "system size x should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(20, se.system_size_y, "system size y should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(20, se.system_size_z, "system size z should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se.num_elements, "num elements should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, se.num_nn_levels, "num nn levels should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(6, se.num_bond_types, "num bond types should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(8, se.num_neighbor_types, "num neighbor types should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12, se.num_transition_vectors,
                                  "num transition vectors should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.0, se.primitive_basis[0][0],
                                     "primitive basis [0][0] should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.0, se.primitive_basis[1][1],
                                     "primitive basis [1][1] should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.0, se.primitive_basis[2][2],
                                     "primitive basis [2][2] should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.01, se.overpotential_ramp_rate,
                                     "overpotential ramp rate should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.0, se.max_overpotential, "max overpotential should restore");

    remove(checkpoint_path);
}

void test_checkpoint_save_and_load_state_and_env(void)
{
    const char checkpoint_path[] = "build/test/checkpoint_state_and_env.bin";
    struct SimulationState ss = {0};
    struct SimulationEnv se = {0};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    // populate both state and env for a combined checkpoint.
    ss.iter = 100ul;
    ss.temperature = 300.0;
    ss.overpotential = 0.5;
    se.flavor = FLAVOR_MC;
    se.rand_seed = 54321u;
    se.system_size_x = 30;
    se.num_elements = 3;

    save_status = checkpoint_save(checkpoint_path, &ss, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "combined state+env checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));
    memset(&se, 0, sizeof(se));

    load_status = checkpoint_load(checkpoint_path, &ss, &se);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "combined state+env checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(100u, ss.iter, "iter should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(300.0, ss.temperature,
                                     "temperature should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_MC, se.flavor,
                                   "flavor should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(30, se.system_size_x,
                                  "system size should restore from combined checkpoint");

    remove(checkpoint_path);
}

void test_checkpoint_save_and_load_atom_round_trip(void)
{
    const char checkpoint_path[] = "build/test/checkpoint_atom_round_trip.bin";
    struct SimulationState ss = {0};
    Atom original_atom = {0};
    Atom *atom_refs[1] = {&original_atom};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    original_atom.type = 7u;
    original_atom.energy = -12.5;
    original_atom.lattice[0] = 4;
    original_atom.lattice[1] = -2;
    original_atom.lattice[2] = 9;
    original_atom.bsradius = 0.42;

    ss.atom_cnt = 1;
    ss.atom_arr = atom_refs;

    save_status = checkpoint_save(checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "atom checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));

    load_status = checkpoint_load(checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "atom checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(1u, (unsigned int)ss.atom_cnt, "atom count should restore");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss.atom_arr, "atom array should be allocated on restore");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&original_atom, ss.atom_arr[0], sizeof(Atom),
                                     "atom should round-trip byte-for-byte for serialized fields");

    free(ss.atom_arr[0]);
    free(ss.atom_arr);
    remove(checkpoint_path);
}

void test_checkpoint_load_missing_file_returns_error(void)
{
    const char checkpoint_path[] = "build/test/does_not_exist_checkpoint.bin";
    CheckpointStatus load_status;

    remove(checkpoint_path);

    load_status = checkpoint_load(checkpoint_path, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, load_status,
                                  "loading a missing checkpoint should fail cleanly");
}

void test_checkpoint_load_corrupted_payload_returns_error(void)
{
    const char checkpoint_path[] = "build/test/checkpoint_corrupted_payload.bin";
    struct SimulationState ss = {0};
    FILE *file;
    int corrupted_byte;
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    ss.iter = 11ul;
    ss.temperature = 273.15;
    ss.overpotential = 0.25;

    save_status = checkpoint_save(checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "baseline checkpoint should save successfully");

    file = fopen(checkpoint_path, "rb+");
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

    load_status = checkpoint_load(checkpoint_path, &ss, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, load_status,
                                  "loading a corrupted checkpoint should fail cleanly");

    remove(checkpoint_path);
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

    UNITY_END();
    return 0;
}
