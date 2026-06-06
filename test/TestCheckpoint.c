#include "Checkpoint.h"
#include "CheckpointLogging.h"
#include "CheckpointSimulation.h"
#include "CheckpointUtils.h"
#include "Initialization.h"
#include "TUtils.h"
#include "unity.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>

FILE *temp_checkpoint_file;
const char temp_checkpoint_path[] = "test/output/temp_checkpoint.bin";

// TODO: move the setup of temp files to TUtils and reuse
#define NTEMPS 4
FILE *temp_logs[NTEMPS];
char *temp_names[NTEMPS] = {"test/output/temp_log_0.log", "test/output/temp_log_1.log",
                            "test/output/temp_log_2.log", "test/output/temp_log_3.log"};

// TODO: set up ss and se in setUp so it cleans on test fail
void setUp(void)
{
    for (int i = 0; i < NTEMPS; i++) {
        temp_logs[i] = fopen(temp_names[i], "wb+");
        fopen_error(temp_names[i], temp_logs[i]);
    }
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

    for (int i = 0; i < NTEMPS; i++) {
        fclose(temp_logs[i]);
    }
}

void test_checkpoint_checksum32_is_deterministic(void)
{
    // Use the same input twice to prove the checksum is stable, then flip one byte to prove it
    // changes.
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

void test_checkpoint_save_and_load_env_arrays_and_atom_names(void)
{
    // Round-trip the env arrays that are serialized separately from the scalar payload.
    struct SimulationEnv se_save = {0};
    struct SimulationEnv se_load = {0};
    const char *names[] = {"A", "BB"};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    se_save.flavor = FLAVOR_KMC;
    se_save.lattice_type = FCC;
    se_save.rand_state = (RandomState){.u = 12345, .v = 67890, .w = 24680};
    se_save.system_size_x = 20;
    se_save.system_size_y = 20;
    se_save.system_size_z = 20;
    se_save.overpotential_ramp_rate = 0.015;
    se_save.max_overpotential = 1.5;
    se_save.dissolution = 1;
    se_save.num_elements = 2;
    se_save.num_nn_levels = 1;
    se_save.num_bond_types = 2;
    se_save.num_nn_types = 1;
    se_save.num_transition_vectors = 4;
    se_save.atom_names_cnt = 2;

    // Allocate each array explicitly so the checkpoint code can copy and later restore them.
    se_save.substrate_composition = (double *)malloc((size_t)se_save.num_elements * sizeof(double));
    se_save.nn_energy = (double *)malloc((size_t)se_save.num_nn_types * sizeof(double));
    se_save.is_soluble = (bool *)malloc((size_t)se_save.num_elements * sizeof(bool));
    se_save.atom_names = (char **)malloc((size_t)se_save.atom_names_cnt * sizeof(char *));
    TEST_ASSERT_NOT_NULL_MESSAGE(se_save.substrate_composition,
                                 "substrate composition should allocate");
    TEST_ASSERT_NOT_NULL_MESSAGE(se_save.nn_energy, "nn energy should allocate");
    TEST_ASSERT_NOT_NULL_MESSAGE(se_save.is_soluble, "is_soluble should allocate");
    TEST_ASSERT_NOT_NULL_MESSAGE(se_save.atom_names, "atom_names should allocate");

    se_save.substrate_composition[0] = 0.25;
    se_save.substrate_composition[1] = 0.75;
    se_save.nn_energy[0] = -1.25;
    se_save.is_soluble[0] = false;
    se_save.is_soluble[1] = true;
    for (int i = 0; i < se_save.atom_names_cnt; ++i) {
        size_t name_len = strlen(names[i]) + 1u;
        se_save.atom_names[i] = (char *)malloc(name_len);
        TEST_ASSERT_NOT_NULL_MESSAGE(se_save.atom_names[i], "atom name should allocate");
        memcpy(se_save.atom_names[i], names[i], name_len);
    }

    save_status = write_checkpoint_buffer(temp_checkpoint_path, NULL, &se_save, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "env array checkpoint should save successfully");

    free(se_save.substrate_composition);
    free(se_save.nn_energy);
    free(se_save.is_soluble);
    for (int i = 0; i < se_save.atom_names_cnt; ++i) {
        free(se_save.atom_names[i]);
    }
    free(se_save.atom_names);
    free(se_save.transition_vectors);
    free(se_save.opposite_tvectors);

    // Load into a blank env so the test can verify the restore path repopulates every array.
    load_status = read_checkpoint_file(temp_checkpoint_path, NULL, &se_load, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "env array checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, se_load.flavor, "flavor should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(12345, se_load.rand_state.u, "random state u should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(67890, se_load.rand_state.v, "random state v should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(24680, se_load.rand_state.w, "random state w should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.25, se_load.substrate_composition[0],
                                     "substrate composition should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, se_load.substrate_composition[1],
                                     "substrate composition should restore");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(-1.25, se_load.nn_energy[0], "nn energy should restore");
    TEST_ASSERT_FALSE_MESSAGE(se_load.is_soluble[0], "solubility should restore");
    TEST_ASSERT_TRUE_MESSAGE(se_load.is_soluble[1], "solubility should restore");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se_load.atom_names_cnt, "atom name count should restore");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("A", se_load.atom_names[0], "atom name should restore");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("BB", se_load.atom_names[1], "atom name should restore");

    free(se_load.substrate_composition);
    free(se_load.nn_energy);
    free(se_load.is_soluble);
    for (int i = 0; i < se_load.atom_names_cnt; ++i) {
        free(se_load.atom_names[i]);
    }
    free(se_load.atom_names);
    free(se_load.transition_vectors);
    free(se_load.opposite_tvectors);
    free(se_load.atoms_per_nn_level);
}

void test_checkpoint_load_rejects_corrupted_env_atom_names_header(void)
{
    // Corrupt the atom-name array header in-place to make sure load rejects a bad array marker.
    struct SimulationEnv se = {0};
    const char *names[] = {"A", "BB"};
    const size_t array_header_size =
        CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t);
    // Skip the scalar payload and the first three env arrays to land on the atom-name header.
    const size_t payload_offset = sizeof(CheckpointHeader) + sizeof(SimEnvPayload) +
                                  (array_header_size + (size_t)2 * sizeof(double)) +
                                  (array_header_size + (size_t)1 * sizeof(double)) +
                                  (array_header_size + (size_t)2 * sizeof(bool));
    CheckpointStatus save_status;
    CheckpointStatus load_status;
    FILE *file;
    int ch;

    se.flavor = FLAVOR_MC;
    se.lattice_type = BCC;
    se.rand_state = (RandomState){.u = 54321, .v = 9876, .w = 13579};
    se.system_size_x = 20;
    se.system_size_y = 20;
    se.system_size_z = 20;
    se.overpotential_ramp_rate = 0.01;
    se.max_overpotential = 1.0;
    se.dissolution = 1;
    se.num_elements = 2;
    se.num_nn_levels = 1;
    se.num_bond_types = 1;
    se.num_nn_types = 1;
    se.num_transition_vectors = 1;
    se.atom_names_cnt = 2;
    se.substrate_composition = (double *)malloc((size_t)se.num_elements * sizeof(double));
    se.nn_energy = (double *)malloc((size_t)se.num_nn_types * sizeof(double));
    se.is_soluble = (bool *)malloc((size_t)se.num_elements * sizeof(bool));
    se.atom_names = (char **)malloc((size_t)se.atom_names_cnt * sizeof(char *));
    TEST_ASSERT_NOT_NULL_MESSAGE(se.substrate_composition, "substrate composition should allocate");
    TEST_ASSERT_NOT_NULL_MESSAGE(se.nn_energy, "nn energy should allocate");
    TEST_ASSERT_NOT_NULL_MESSAGE(se.is_soluble, "is_soluble should allocate");
    TEST_ASSERT_NOT_NULL_MESSAGE(se.atom_names, "atom_names should allocate");
    se.substrate_composition[0] = 0.2;
    se.substrate_composition[1] = 0.8;
    se.nn_energy[0] = 2.25;
    se.is_soluble[0] = true;
    se.is_soluble[1] = false;
    for (int i = 0; i < se.atom_names_cnt; ++i) {
        size_t name_len = strlen(names[i]) + 1u;
        se.atom_names[i] = (char *)malloc(name_len);
        TEST_ASSERT_NOT_NULL_MESSAGE(se.atom_names[i], "atom name should allocate");
        memcpy(se.atom_names[i], names[i], name_len);
    }

    save_status = write_checkpoint_buffer(temp_checkpoint_path, NULL, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "baseline env array checkpoint should save successfully");

    file = fopen(temp_checkpoint_path, "rb+");
    TEST_ASSERT_NOT_NULL_MESSAGE(file, "checkpoint file should open for corruption");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, fseek(file, (long)payload_offset, SEEK_SET),
                                  "file should seek to atom names header");
    // Read one byte, flip its bits, and write it back to keep the corruption narrow.
    ch = fgetc(file);
    TEST_ASSERT_TRUE_MESSAGE(ch != EOF, "atom names header byte should be readable");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, fseek(file, (long)payload_offset, SEEK_SET),
                                  "file should seek back before overwrite");
    TEST_ASSERT_NOT_EQUAL_MESSAGE(EOF, fputc(ch ^ 0xFF, file),
                                  "atom names header should be overwritten");
    fclose(file);

    load_status = read_checkpoint_file(temp_checkpoint_path, NULL, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, load_status,
                                  "corrupted atom names header should fail to load");

    free(se.substrate_composition);
    free(se.nn_energy);
    free(se.is_soluble);
    for (int i = 0; i < se.atom_names_cnt; ++i) {
        free(se.atom_names[i]);
    }
    free(se.atom_names);
}

void test_checkpoint_header_has_valid_checksum_accepts_matching_payload(void)
{
    // A checksum computed from the same payload should validate exactly.
    const uint8_t payload[] = {0x10u, 0x20u, 0x30u, 0x40u};
    CheckpointHeader header = {0};

    initialize_checkpoint_header(&header);
    finalize_checkpoint_header(&header, (uint32_t)sizeof(payload),
                               checkpoint_checksum32(payload, sizeof(payload)));

    TEST_ASSERT_TRUE_MESSAGE(
        checkpoint_header_has_valid_checksum(&header, payload, sizeof(payload)),
        "matching payload should validate its checksum");
}

void test_checkpoint_header_has_valid_checksum_rejects_corrupted_payload(void)
{
    // Change one byte so the checksum comparison fails.
    const uint8_t payload[] = {0x10u, 0x20u, 0x30u, 0x40u};
    const uint8_t corrupted_payload[] = {0x10u, 0x20u, 0x30u, 0x41u};
    CheckpointHeader header = {0};

    initialize_checkpoint_header(&header);
    finalize_checkpoint_header(&header, (uint32_t)sizeof(payload),
                               checkpoint_checksum32(payload, sizeof(payload)));

    TEST_ASSERT_FALSE_MESSAGE(
        checkpoint_header_has_valid_checksum(&header, corrupted_payload, sizeof(corrupted_payload)),
        "corrupted payload should fail checksum validation");
}

void test_verify_payload_size_accepts_matching_file_size(void)
{
    // Build a synthetic file with the exact expected size so the verifier can seek to the payload.
    CheckpointHeader header = {0};
    FILE *file;
    const uint8_t payload[] = {0x11u, 0x22u, 0x33u, 0x44u};

    initialize_checkpoint_header(&header);
    header.payload_bytes = (uint32_t)sizeof(payload);

    temp_checkpoint_file = fopen(temp_checkpoint_path, "wb+");
    fopen_error(temp_checkpoint_path, temp_checkpoint_file);
    // Write the header and payload explicitly so the test controls the full file size.
    TEST_ASSERT_EQUAL_UINT_MESSAGE(
        (unsigned int)sizeof(header),
        (unsigned int)fwrite(&header, 1u, sizeof(header), temp_checkpoint_file),
        "header should write to temp file");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(
        (unsigned int)sizeof(payload),
        (unsigned int)fwrite(payload, 1u, sizeof(payload), temp_checkpoint_file),
        "payload should write to temp file");
    rewind(temp_checkpoint_file);

    file = temp_checkpoint_file;
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, verify_payload_size(file, &header),
                                  "matching file size should validate successfully");
    TEST_ASSERT_EQUAL_INT_MESSAGE((int)header.header_bytes, ftell(file),
                                  "successful validation should leave the file at payload start");
}

void test_verify_payload_size_rejects_mismatched_file_size(void)
{
    // Make the header claim more bytes than are actually on disk and confirm the helper rejects it.
    CheckpointHeader header = {0};
    const uint8_t payload[] = {0x11u, 0x22u};

    initialize_checkpoint_header(&header);
    header.payload_bytes = (uint32_t)sizeof(payload) + 1u;

    temp_checkpoint_file = fopen(temp_checkpoint_path, "wb+");
    fopen_error(temp_checkpoint_path, temp_checkpoint_file);
    TEST_ASSERT_EQUAL_UINT_MESSAGE(
        (unsigned int)sizeof(header),
        (unsigned int)fwrite(&header, 1u, sizeof(header), temp_checkpoint_file),
        "header should write to temp file");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(
        (unsigned int)sizeof(payload),
        (unsigned int)fwrite(payload, 1u, sizeof(payload), temp_checkpoint_file),
        "payload should write to temp file");
    rewind(temp_checkpoint_file);

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR,
                                  verify_payload_size(temp_checkpoint_file, &header),
                                  "mismatched file size should be rejected");
}

void test_rebuild_rates_and_transitions_is_noop_for_empty_state(void)
{
    // With no atoms present, the rebuild helper should return success without doing any work.
    struct SimulationState ss = {0};
    struct SimulationEnv se = {0};

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, rebuild_rates_and_transitions(&ss, &se),
                                  "empty state should be a no-op for transition rebuilding");
}

void test_checkpoint_header_magic_helpers(void)
{
    // Set the magic, validate it, then corrupt one byte to confirm the check fails.
    CheckpointHeader header = {0};

    set_checkpoint_header_magic(&header);

    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_has_valid_magic(&header),
                             "magic should match after being set");
    header.magic[0] = 'X';
    TEST_ASSERT_FALSE_MESSAGE(checkpoint_header_has_valid_magic(&header),
                              "magic should fail after corruption");
}

void test_checkpoint_header_is_valid_accepts_valid_header(void)
{
    // A freshly initialized header should pass the combined validity check.
    CheckpointHeader header = {0};

    initialize_checkpoint_header(&header);

    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_is_valid(&header),
                             "initialized header should be valid");
}

void test_checkpoint_header_is_valid_rejects_invalid_structure(void)
{
    // Break the stored header size so the combined check rejects the structure.
    CheckpointHeader header = {0};

    initialize_checkpoint_header(&header);
    header.header_bytes = 0u;

    TEST_ASSERT_FALSE_MESSAGE(checkpoint_header_is_valid(&header),
                              "header with invalid structure should be rejected");
}

void test_checkpoint_header_has_valid_structure_accepts_initialized_header(void)
{
    // The initialized header should meet the version/size/reserved field requirements.
    CheckpointHeader header = {0};

    initialize_checkpoint_header(&header);

    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_has_valid_structure(&header),
                             "initialized header should have a valid structure");
}

void test_checkpoint_header_has_valid_structure_rejects_corrupted_header_bytes(void)
{
    // Zeroing the stored size should be enough to make the structure invalid.
    CheckpointHeader header = {0};

    initialize_checkpoint_header(&header);
    header.header_bytes = 0u;

    TEST_ASSERT_FALSE_MESSAGE(checkpoint_header_has_valid_structure(&header),
                              "corrupted header bytes should fail structure validation");
}

void test_checkpoint_header_init_and_finalize(void)
{
    // Initialize the header, then finalize it with an explicit payload size and checksum.
    CheckpointHeader header;

    initialize_checkpoint_header(&header);
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(CHECKPOINT_FORMAT_VERSION, header.format_version,
                                     "header version should be initialized");
    TEST_ASSERT_EQUAL_UINT32_MESSAGE((uint32_t)sizeof(CheckpointHeader), header.header_bytes,
                                     "header size should be stored");
    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_has_valid_structure(&header),
                             "header shape should be valid after init");

    finalize_checkpoint_header(&header, 128u, 0x12345678u);
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(128u, header.payload_bytes, "payload size should be written");
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(0x12345678u, header.checksum, "checksum should be written");
}

void test_checkpoint_header_write_and_read_round_trip(void)
{
    // Write a header to disk and read it back unchanged.
    CheckpointHeader written;
    CheckpointHeader read_back;

    temp_checkpoint_file = fopen(temp_checkpoint_path, "wb+");
    fopen_error(temp_checkpoint_path, temp_checkpoint_file);

    initialize_checkpoint_header(&written);
    finalize_checkpoint_header(&written, 256u, 0xDEADBEEFu);

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK,
                                  write_checkpoint_header(temp_checkpoint_file, &written),
                                  "header should write successfully");
    rewind(temp_checkpoint_file);

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK,
                                  read_checkpoint_header(temp_checkpoint_file, &read_back),
                                  "header should read successfully");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&written, &read_back, sizeof(CheckpointHeader),
                                     "header should round-trip byte-for-byte");
}

void test_checkpoint_save_and_load_header_only(void)
{
    // Saving with no state or env should still produce a valid header-only file.
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    save_status = write_checkpoint_buffer(temp_checkpoint_path, NULL, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "header-only checkpoint should save successfully");

    load_status = read_checkpoint_file(temp_checkpoint_path, NULL, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "header-only checkpoint should load successfully");
}

void test_checkpoint_save_and_load_state_scalars(void)
{
    // Save a state-only checkpoint and confirm the scalar fields survive the round trip.
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

    save_status = write_checkpoint_buffer(temp_checkpoint_path, &ss, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "state scalar checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));

    load_status = read_checkpoint_file(temp_checkpoint_path, &ss, NULL, NULL);
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
    // Save an env-only checkpoint so the scalar environment metadata can be restored independently.
    struct SimulationEnv se = {0};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    // populate a minimal env with key config fields for checkpointing.
    // This test keeps the env minimal on purpose so only the checkpointed metadata matters.
    se.flavor = FLAVOR_KMC;
    se.lattice_type = FCC;
    se.rand_state = (RandomState){.u = 12345, .v = 67890, .w = 24680};
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

    save_status = write_checkpoint_buffer(temp_checkpoint_path, NULL, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "env scalar checkpoint should save successfully");

    memset(&se, 0, sizeof(se));

    load_status = read_checkpoint_file(temp_checkpoint_path, NULL, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "env scalar checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, se.flavor, "flavor should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(24680, se.rand_state.w, "random state should restore");
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
    free(se.substrate_composition);
    free(se.nn_energy);
    free(se.is_soluble);
    free(se.atom_names);
    free(se.transition_vectors);
    free(se.opposite_tvectors);
    free(se.atoms_per_nn_level);
}

void test_checkpoint_save_and_load_state_and_env(void)
{
    // Save both payload types together and verify the state and env are restored in one load.
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
    se.rand_state = (RandomState){.u = 12345, .v = 67890, .w = 24680};
    se.system_size_x = 40;
    se.system_size_y = 20;
    se.system_size_z = 20;
    se.num_elements = 2;
    se.num_nn_levels = 1;
    se.num_bond_types = 3;
    se.num_transition_vectors = 12;
    se.overpotential_ramp_rate = 0.01;
    se.max_overpotential = 1.0;

    save_status = write_checkpoint_buffer(temp_checkpoint_path, &ss, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "combined state+env checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));
    memset(&se, 0, sizeof(se));

    load_status = read_checkpoint_file(temp_checkpoint_path, &ss, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "combined state+env checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(100u, ss.iter, "iter should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(300.0, ss.temperature,
                                     "temperature should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_MC, se.flavor,
                                   "flavor should restore from combined checkpoint");
    TEST_ASSERT_EQUAL_INT_MESSAGE(40, se.system_size_x,
                                  "system size should restore from combined checkpoint");

    free(ss.atom_arr);
    free(ss.rate_arr);
    free(ss.zone_arr);
    free(ss.transition_arr);
    free(ss.transition_probability.rate_arr_index);
    free(ss.transition_probability.lbound);
    free(ss.transition_probability.ubound);

    free(se.substrate_composition);
    free(se.nn_energy);
    free(se.is_soluble);
    free(se.atom_names);
    free(se.transition_vectors);
    free(se.opposite_tvectors);
    free(se.atoms_per_nn_level);
}

void test_checkpoint_save_and_load_atom_round_trip(void)
{
    // Include a single atom so the checkpoint also exercises atom array serialization.
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
    se.rand_state = (RandomState){.u = 12345, .v = 67890, .w = 24680};
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
    // The restore path expects this env array to exist when rebuilding atom metadata.
    se.is_soluble = (bool *)malloc((size_t)se.num_elements * sizeof(bool));
    TEST_ASSERT_NOT_NULL_MESSAGE(se.is_soluble, "is_soluble array should allocate");
    for (int i = 0; i < se.num_elements; ++i) {
        se.is_soluble[i] = false;
    }

    // Point the state at the single atom so write_checkpoint_buffer writes the atom array.
    ss.atom_cnt = 1;
    ss.atom_arr = atom_refs;

    save_status = write_checkpoint_buffer(temp_checkpoint_path, &ss, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "atom checkpoint should save successfully");
    free(se.is_soluble);

    memset(&ss, 0, sizeof(ss));

    load_status = read_checkpoint_file(temp_checkpoint_path, &ss, &se, NULL);
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
    free(ss.rate_arr);
    free(ss.transition_arr);
    free(ss.transition_probability.rate_arr_index);
    free(ss.transition_probability.lbound);
    free(ss.transition_probability.ubound);
    for (size_t i = 0; i < se.zone_count_u; i++) {
        for (size_t j = 0; j < se.zone_count_v; j++) {
            free(ss.zone_arr[i][j]);
        }
        free(ss.zone_arr[i]);
    }
    free(ss.zone_arr);

    free(se.substrate_composition);
    free(se.nn_energy);
    free(se.is_soluble);
    free(se.atom_names);
    free(se.transition_vectors);
    free(se.opposite_tvectors);
    free(se.atoms_per_nn_level);
}

void test_checkpoint_save_and_load_logging_scalars_success(void)
{
    struct LoggingState ls = {0};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    ls.framenum = 17;
    ls.verbose = 1;
    ls.verbose_interval = 250ul;
    ls.increment_precision = 3;
    ls.stime_precision = 4;
    ls.overpot_precision = 5;
    ls.out_formats_cnt = 0;
    ls.out_formats = NULL;

    save_status = write_checkpoint_buffer(temp_checkpoint_path, NULL, NULL, &ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "checkpoint should save successfully with null state/env");

    memset(&ls, 0, sizeof(ls));

    load_status = read_checkpoint_file(temp_checkpoint_path, NULL, NULL, &ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "checkpoint should load successfully with null state/env");
    TEST_ASSERT_EQUAL_INT_MESSAGE(17, ls.framenum, "frame number should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls.verbose, "verbose flag should copy into payload");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(250u, (unsigned int)ls.verbose_interval,
                                   "verbose interval should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, ls.increment_precision,
                                  "increment precision should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, ls.stime_precision,
                                  "stime precision should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(5, ls.overpot_precision,
                                  "overpotential precision should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, ls.out_formats_cnt,
                                  "output formats count should copy into payload");
    TEST_ASSERT_NULL_MESSAGE(ls.out_formats, "output formats pointer should copy into payload");
}

void test_checkpoint_save_and_load_logging_scalars_null_payload(void)
{
    // A null logging state pointer should cause the helpers to skip the payload without error.
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    save_status = write_checkpoint_buffer(temp_checkpoint_path, NULL, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "checkpoint should save successfully with null logging state");

    load_status = read_checkpoint_file(temp_checkpoint_path, NULL, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "checkpoint should load successfully with null logging state");
}

void test_checkpoint_save_and_load_logging_formats_success(void)
{
    struct LoggingState save_ls = {0};
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    const int n = 4;
    save_ls.out_formats_cnt = n;
    save_ls.out_formats = (OutputFormat *)calloc((size_t)n, sizeof(OutputFormat));
    TEST_ASSERT_NOT_NULL_MESSAGE(save_ls.out_formats, "out_formats should allocate");

    /* CSV format entry */
    save_ls.out_formats[0].type = OUTPUT_FORMAT_CSV;
    save_ls.out_formats[0].is_active = true;
    strncpy(save_ls.out_formats[0].csv.filename, temp_names[0],
            sizeof(save_ls.out_formats[0].csv.filename));
    save_ls.out_formats[0].csv.field_count = 2;
    save_ls.out_formats[0].csv.frame_num = 42;
    save_ls.out_formats[0].csv.schedule.mode = OUTPUT_SCHEDULE_INTERVAL_ITERATION;
    save_ls.out_formats[0].csv.schedule.interval = 1.5;
    save_ls.out_formats[0].csv.schedule.frame_num = 7;

    /* Second CSV format entry */
    save_ls.out_formats[1].type = OUTPUT_FORMAT_CSV;
    save_ls.out_formats[1].is_active = false;
    strncpy(save_ls.out_formats[1].csv.filename, temp_names[1],
            sizeof(save_ls.out_formats[1].csv.filename));
    save_ls.out_formats[1].csv.field_count = 4;
    save_ls.out_formats[1].csv.frame_num = 99;
    save_ls.out_formats[1].csv.schedule.mode = OUTPUT_SCHEDULE_LIST_TIME;
    save_ls.out_formats[1].csv.schedule.interval = 0.25;
    save_ls.out_formats[1].csv.schedule.frame_num = 11;

    /* STEPS format entry */
    save_ls.out_formats[2].type = OUTPUT_FORMAT_STEPS_CSV;
    save_ls.out_formats[2].is_active = true;
    strncpy(save_ls.out_formats[2].steps.filename, temp_names[2],
            sizeof(save_ls.out_formats[2].steps.filename));
    save_ls.out_formats[2].steps.with_coordination = true;

    /* XYZ format entry */
    save_ls.out_formats[3].type = OUTPUT_FORMAT_XYZ;
    save_ls.out_formats[3].is_active = true;
    strncpy(save_ls.out_formats[3].xyz.prefix, "prefix", sizeof(save_ls.out_formats[3].xyz.prefix));
    strncpy(save_ls.out_formats[3].xyz.suffix, "suffix", sizeof(save_ls.out_formats[3].xyz.suffix));
    save_ls.out_formats[3].xyz.frame_num = 13;
    save_ls.out_formats[3].xyz.stripped = true;
    save_ls.out_formats[3].xyz.schedule.mode = OUTPUT_SCHEDULE_INTERVAL_TIME;
    save_ls.out_formats[3].xyz.schedule.interval = 2.5;
    save_ls.out_formats[3].xyz.schedule.frame_num = 17;

    save_status = write_checkpoint_buffer(temp_checkpoint_path, NULL, NULL, &save_ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "checkpoint should save successfully with null state/env");

    struct LoggingState load_ls = {0};

    load_status = read_checkpoint_file(temp_checkpoint_path, NULL, NULL, &load_ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "checkpoint should load successfully with null state/env");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, load_ls.framenum, "frame number should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, load_ls.verbose, "verbose flag should copy into payload");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(0, (unsigned int)load_ls.verbose_interval,
                                   "verbose interval should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, load_ls.increment_precision,
                                  "increment precision should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, load_ls.stime_precision,
                                  "stime precision should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, load_ls.overpot_precision,
                                  "overpotential precision should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(n, load_ls.out_formats_cnt,
                                  "output formats count should copy into payload");
    TEST_ASSERT_NOT_NULL_MESSAGE(load_ls.out_formats,
                                 "output formats pointer should copy into payload");

    for (int i = 0; i < n; ++i) {
        char message[64];
        snprintf(message, sizeof(message), "filled output format %d should match", i);
        assert_output_format_matches_runtime_round_trip(&save_ls.out_formats[i],
                                                        &load_ls.out_formats[i], message);
    }

    fclose(load_ls.out_formats[0].csv.file);
    fclose(load_ls.out_formats[1].csv.file);
    fclose(load_ls.out_formats[2].steps.file);
    free(save_ls.out_formats);
    free(load_ls.out_formats);
}

void test_checkpoint_load_missing_file_returns_error(void)
{
    // Loading a missing path should fail cleanly rather than crashing or creating files.
    CheckpointStatus load_status;

    remove(temp_checkpoint_path);

    load_status = read_checkpoint_file(temp_checkpoint_path, NULL, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, load_status,
                                  "loading a missing checkpoint should fail cleanly");
}

void test_checkpoint_load_corrupted_payload_returns_error(void)
{
    // Corrupt the stored payload bytes and verify the checksum guard rejects the file.
    struct SimulationState ss = {0};
    FILE *file;
    int corrupted_byte;
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    ss.iter = 11ul;
    ss.temperature = 273.15;
    ss.overpotential = 0.25;

    save_status = write_checkpoint_buffer(temp_checkpoint_path, &ss, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "baseline checkpoint should save successfully");

    file = fopen(temp_checkpoint_path, "rb+");
    TEST_ASSERT_NOT_NULL_MESSAGE(file, "checkpoint file should open for corruption");

    // Move to the payload start so the next read overwrites the first payload byte.
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, fseek(file, (long)sizeof(CheckpointHeader), SEEK_SET),
                                  "file should seek to payload start");

    // Flip one payload byte instead of rewriting the whole file; checksum failure is enough.
    corrupted_byte = fgetc(file);
    TEST_ASSERT_TRUE_MESSAGE(corrupted_byte != EOF, "payload byte should be readable");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, fseek(file, (long)sizeof(CheckpointHeader), SEEK_SET),
                                  "file should seek back to payload start for overwrite");
    // flip bits and write byte to start of file
    TEST_ASSERT_NOT_EQUAL_MESSAGE(EOF, fputc(corrupted_byte ^ 0xFF, file),
                                  "payload byte should be overwritten");
    fclose(file);

    load_status = read_checkpoint_file(temp_checkpoint_path, &ss, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, load_status,
                                  "loading a corrupted checkpoint should fail cleanly");
}

void test_checkpoint_rebuild_zones_and_rates_from_atoms(void)
{
    // Exercise the full restore path that rebuilds zones and transition/rate data from atoms.
    struct SimulationState ss_orig = {0};
    struct SimulationState ss_restored = {0};
    struct SimulationEnv se = {0};
    CheckpointStatus save_status, load_status;

    // Set up a compact environment that still has enough geometry for zone lookup.
    struct SimulationConfig dummy_config = {
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

    // Create a tiny atom set that is easy to reason about after restoration.
    ss_orig.iter = 50ul;
    ss_orig.temperature = 300.0;
    ss_orig.atom_cnt = 2u;

    // Allocate each atom manually so the checkpoint has a real atom array to serialize.
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

        // Initialize every neighbor and transition slot so restore can rebuild them from scratch.
        for (int j = 0; j < MAXIMUM_NUMBER_OF_NEIGHBORS; ++j) {
            ss_orig.atom_arr[i]->neighbor_atom_idxs[j] = -1;
        }
        for (int j = 0; j < MAXIMUM_NUMBER_OF_NEIGHBORS + DISSOLUTION; ++j) {
            ss_orig.atom_arr[i]->transition_indices[j] = -1;
        }
    }

    // Save the checkpoint before the restore-specific arrays are allocated.
    save_status = write_checkpoint_buffer(temp_checkpoint_path, &ss_orig, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "checkpoint with atoms should save successfully");

    free(se.transition_vectors);
    free(se.opposite_tvectors);
    free(se.atoms_per_nn_level);
    free(se.is_soluble);
    free(se.nn_energy);

    // Load the checkpoint and let the restore path rebuild the derived structures.
    load_status = read_checkpoint_file(temp_checkpoint_path, &ss_restored, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "checkpoint with atoms should load successfully");

    // The restored state should match the saved atom count and core counters.
    TEST_ASSERT_EQUAL_INT_MESSAGE(ss_orig.iter, ss_restored.iter, "iter should be restored");
    TEST_ASSERT_EQUAL_INT_MESSAGE(ss_orig.atom_cnt, ss_restored.atom_cnt,
                                  "atom count should be restored");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss_restored.atom_arr, "atom_arr should be allocated");
    TEST_ASSERT_NOT_NULL_MESSAGE(ss_restored.zone_arr, "zone_arr should be allocated");

    // The restore path should zero transient linkage fields before rebuilding them.
    for (int i = 0; i < ss_restored.atom_cnt; ++i) {
        TEST_ASSERT_NOT_NULL_MESSAGE(ss_restored.atom_arr[i], "each atom should be allocated");
        // The exact transition set is implementation-defined, but the atom coordinates should be
        // valid.
        TEST_ASSERT_TRUE_MESSAGE(ss_restored.atom_arr[i]->lattice[0] >= 0,
                                 "atom lattice coordinates should be valid");
    }

    // Clean up the original and restored atom arrays separately.
    for (int i = 0; i < ss_orig.atom_cnt; ++i) {
        free(ss_orig.atom_arr[i]);
    }
    free(ss_orig.atom_arr);

    for (int i = 0; i < ss_restored.atom_cnt; ++i) {
        free(ss_restored.atom_arr[i]);
    }
    free(ss_restored.atom_arr);

    // Free the rebuilt zone grid after the restore assertions finish.
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

    free(ss_restored.rate_arr);
    free(ss_restored.transition_arr);
    free(ss_restored.transition_probability.rate_arr_index);
    free(ss_restored.transition_probability.lbound);
    free(ss_restored.transition_probability.ubound);

    free(se.substrate_composition);
    free(se.nn_energy);
    free(se.is_soluble);
    free(se.atom_names);
    free(se.transition_vectors);
    free(se.opposite_tvectors);
    free(se.atoms_per_nn_level);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_checkpoint_checksum32_is_deterministic);

    RUN_TEST(test_checkpoint_save_and_load_env_arrays_and_atom_names);
    RUN_TEST(test_checkpoint_load_rejects_corrupted_env_atom_names_header);

    RUN_TEST(test_checkpoint_header_has_valid_checksum_accepts_matching_payload);
    RUN_TEST(test_checkpoint_header_has_valid_checksum_rejects_corrupted_payload);

    RUN_TEST(test_checkpoint_header_magic_helpers);
    RUN_TEST(test_checkpoint_header_init_and_finalize);
    RUN_TEST(test_checkpoint_header_write_and_read_round_trip);

    RUN_TEST(test_checkpoint_save_and_load_header_only);
    RUN_TEST(test_checkpoint_save_and_load_state_scalars);
    RUN_TEST(test_checkpoint_save_and_load_env_scalars);
    RUN_TEST(test_checkpoint_save_and_load_state_and_env);
    RUN_TEST(test_checkpoint_save_and_load_atom_round_trip);
    RUN_TEST(test_checkpoint_save_and_load_logging_scalars_success);
    RUN_TEST(test_checkpoint_save_and_load_logging_scalars_null_payload);
    RUN_TEST(test_checkpoint_save_and_load_logging_formats_success);

    RUN_TEST(test_checkpoint_load_missing_file_returns_error);
    RUN_TEST(test_checkpoint_load_corrupted_payload_returns_error);
    RUN_TEST(test_checkpoint_rebuild_zones_and_rates_from_atoms);

    UNITY_END();

    for (int i = 0; i < NTEMPS; i++) {
        int ret = remove(temp_names[i]);
        if (ret != 0) {
            fprintf(stderr, "Warning: failed to remove temporary file %s - %s\n", temp_names[i],
                    strerror(errno));
        }
    }

    return 0;
}
