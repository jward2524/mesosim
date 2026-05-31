#include "Checkpoint.h"
#include "Initialization.h"
#include "TUtils.h"
#include "unity.h"
#include <stdlib.h>
#include <string.h>

FILE *temp_checkpoint_file;
const char temp_checkpoint_path[] = "build/test/temp_checkpoint.bin";

// TODO: set up ss and se in setUp so it cleans on test fail
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

static void build_array_header(uint8_t *buffer, uint16_t flag, uint32_t n)
{
    memcpy(buffer, &((const uint16_t){CHECKPOINT_ARRAY_MAGIC}), sizeof(uint16_t));
    memcpy(buffer + sizeof(uint16_t), &flag, sizeof(uint16_t));
    memcpy(buffer + sizeof(uint16_t) + sizeof(uint16_t), &n, sizeof(uint32_t));
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

void test_fill_state_payload_copies_selected_state_scalars(void)
{
    // Seed the state with values that should survive the copy into the compact payload.
    struct SimulationState ss = {0};
    CheckpointStatePayload payload = {0};

    ss.iter = 17ul;
    ss.mmc_steps = 23ul;
    ss.final_iteration = 99ul;
    ss.run_stime = 4.5;
    ss.simulation_should_kill_itself = true;
    ss.elapsed_stime = 1.25;
    ss.sim_end_type = SIM_END_BY_ITERATIONS;
    ss.frequency_sum = 0.75;
    ss.total_internal_energy = -8.5;
    ss.temperature = 312.0;
    ss.overpotential = 0.42;
    ss.total_atoms_dissolved = 3;

    fill_state_payload(&payload, &ss);

    // Mutate the source after filling the payload to confirm the payload owns its copied values.
    ss.iter = 999ul;

    TEST_ASSERT_EQUAL_UINT_MESSAGE(17u, payload.iter, "iter should copy into the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(ss.iter, payload.iter,
                                       "payload should keep the copied iter value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(23u, payload.mmc_steps,
                                   "mmc steps should copy into the payload");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(99u, payload.final_iteration,
                                   "final iteration should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(4.5, payload.run_stime,
                                     "run time should copy into the payload");
    TEST_ASSERT_TRUE_MESSAGE(payload.simulation_should_kill_itself,
                             "kill-itself flag should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.25, payload.elapsed_stime,
                                     "elapsed time should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, payload.sim_end_type,
                                  "end type should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, payload.frequency_sum,
                                     "frequency sum should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(-8.5, payload.total_internal_energy,
                                     "internal energy should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(312.0, payload.temperature,
                                     "temperature should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.42, payload.overpotential,
                                     "overpotential should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, payload.total_atoms_dissolved,
                                  "dissolved atom count should copy into the payload");
}

void test_apply_state_payload_to_simstate_restores_selected_state_scalars(void)
{
    // Start with a payload that represents the saved checkpoint and a state that should be
    // overwritten.
    CheckpointStatePayload payload = {
        .iter = 11ul,
        .mmc_steps = 22ul,
        .final_iteration = 33ul,
        .run_stime = 6.25,
        .simulation_should_kill_itself = false,
        .elapsed_stime = 2.5,
        .sim_end_type = SIM_END_BY_ITERATIONS,
        .frequency_sum = 0.125,
        .total_internal_energy = -2.75,
        .temperature = 290.0,
        .overpotential = 0.9,
        .total_atoms_dissolved = 6,
    };
    struct SimulationState ss = {
        .iter = 555ul,
        .mmc_steps = 666ul,
        .final_iteration = 777ul,
        .run_stime = 8.0,
        .simulation_should_kill_itself = true,
        .elapsed_stime = 9.0,
        .sim_end_type = -1,
        .frequency_sum = 9.5,
        .total_internal_energy = 12.0,
        .temperature = 333.0,
        .overpotential = 1.5,
        .total_atoms_dissolved = 44,
        .rate_cnt = 101,
        .transition_cnt = 202,
    };

    // Apply the payload back onto the live state and verify only the intended fields change.
    apply_state_payload_to_simstate(&payload, &ss);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(11u, ss.iter, "iter should restore from the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(555u, ss.iter, "iter should overwrite the old value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(22u, ss.mmc_steps, "mmc steps should restore from payload");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(33u, ss.final_iteration,
                                   "final iteration should restore from payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(6.25, ss.run_stime, "run time should restore from payload");
    TEST_ASSERT_FALSE_MESSAGE(ss.simulation_should_kill_itself,
                              "kill-itself flag should restore from payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(2.5, ss.elapsed_stime,
                                     "elapsed time should restore from payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(SIM_END_BY_ITERATIONS, ss.sim_end_type,
                                  "end type should restore from payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.125, ss.frequency_sum,
                                     "frequency sum should restore from payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(-2.75, ss.total_internal_energy,
                                     "internal energy should restore from payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(290.0, ss.temperature,
                                     "temperature should restore from payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.9, ss.overpotential,
                                     "overpotential should restore from payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(6, ss.total_atoms_dissolved,
                                  "dissolved atom count should restore from payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(101, ss.rate_cnt,
                                  "rate count should not be modified by the helper");
    TEST_ASSERT_EQUAL_INT_MESSAGE(202, ss.transition_cnt,
                                  "transition count should not be modified by the helper");
}

void test_fill_env_payload_copies_selected_env_scalars(void)
{
    // This test covers the scalar environment metadata that is stored directly in the payload.
    struct SimulationEnv se = {0};
    CheckpointEnvPayload payload = {0};

    se.flavor = FLAVOR_KMC;
    se.rand_seed = 12345u;
    se.overpotential_ramp_rate = 0.0125;
    se.max_overpotential = 1.25;
    se.system_size_x = 40;
    se.system_size_y = 41;
    se.system_size_z = 42;
    se.num_elements = 3;
    se.num_nn_levels = 2;
    se.num_bond_types = 7;
    se.num_nn_types = 9;
    se.num_transition_vectors = 12;
    se.dissolution = 1;
    se.atom_names_cnt = 4;
    se.lattice_type = FCC;

    fill_env_payload(&payload, &se);

    // Change the source after copying so the payload value is the one under test.
    se.flavor = FLAVOR_MC;

    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, payload.flavor,
                                   "flavor should copy into the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(se.flavor, payload.flavor,
                                       "payload should keep the copied flavor value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(12345u, payload.rand_seed,
                                   "random seed should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.0125, payload.overpotential_ramp_rate,
                                     "ramp rate should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.25, payload.max_overpotential,
                                     "max overpotential should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(40, payload.system_size_x,
                                  "system size x should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(41, payload.system_size_y,
                                  "system size y should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(42, payload.system_size_z,
                                  "system size z should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, payload.num_elements,
                                  "num elements should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, payload.num_nn_levels,
                                  "num nn levels should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(7, payload.num_bond_types,
                                  "num bond types should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(9, payload.num_nn_types,
                                  "num nn types should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12, payload.num_transition_vectors,
                                  "num transition vectors should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, payload.dissolution,
                                  "dissolution flag should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, payload.atom_names_cnt,
                                  "atom name count should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(FCC, payload.lattice_type,
                                  "lattice type should copy into the payload");
}

void test_apply_env_payload_to_config_restores_selected_env_scalars(void)
{
    // The config begins with different values so the restore step has visible work to do.
    CheckpointEnvPayload payload = {
        .flavor = FLAVOR_MC,
        .rand_seed = 54321u,
        .overpotential_ramp_rate = 0.02,
        .max_overpotential = 2.5,
        .system_size_x = 10,
        .system_size_y = 11,
        .system_size_z = 12,
        .num_elements = 5,
        .num_nn_levels = 3,
        .num_bond_types = 8,
        .num_nn_types = 15,
        .num_transition_vectors = 24,
        .dissolution = 0,
        .atom_names_cnt = 2,
        .lattice_type = BCC,
    };
    struct SimulationConfig config = {
        .flavor = FLAVOR_KMC,
        .rand_seed = 1u,
        .overpotential_ramp_rate = 9.0,
        .max_overpotential = 9.5,
        .system_size_x = 99,
        .system_size_y = 98,
        .system_size_z = 97,
        .num_elements = 1,
        .num_nn_levels = 1,
        .num_bond_types = 1,
        .num_nn_types = 1,
        .dissolution = 1,
        .atom_names_cnt = 9,
        .lattice_type = FCC,
    };

    // Apply the saved scalar metadata to the runtime config.
    apply_env_payload_to_config(&payload, &config);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_MC, config.flavor,
                                   "flavor should restore from the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, config.flavor,
                                       "flavor should overwrite the old value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(54321u, config.rand_seed,
                                   "random seed should restore from the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.02, config.overpotential_ramp_rate,
                                     "ramp rate should restore from the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(2.5, config.max_overpotential,
                                     "max overpotential should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(10, config.system_size_x,
                                  "system size x should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(11, config.system_size_y,
                                  "system size y should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12, config.system_size_z,
                                  "system size z should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(5, config.num_elements,
                                  "num elements should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, config.num_nn_levels,
                                  "num nn levels should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(8, config.num_bond_types,
                                  "num bond types should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(15, config.num_nn_types,
                                  "num nn types should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, config.dissolution,
                                  "dissolution should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, config.atom_names_cnt,
                                  "atom name count should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(BCC, config.lattice_type,
                                  "lattice type should restore from the payload");
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
    se_save.rand_seed = 24680u;
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

    save_status = checkpoint_save(temp_checkpoint_path, NULL, &se_save, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "env array checkpoint should save successfully");

    // Load into a blank env so the test can verify the restore path repopulates every array.
    load_status = checkpoint_load(temp_checkpoint_path, NULL, &se_load, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, load_status,
                                  "env array checkpoint should load successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, se_load.flavor, "flavor should restore");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(24680u, se_load.rand_seed, "seed should restore");
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

    free(se_save.substrate_composition);
    free(se_save.nn_energy);
    free(se_save.is_soluble);
    for (int i = 0; i < se_save.atom_names_cnt; ++i) {
        free(se_save.atom_names[i]);
    }
    free(se_save.atom_names);
}

void test_checkpoint_load_rejects_corrupted_env_atom_names_header(void)
{
    // Corrupt the atom-name array header in-place to make sure load rejects a bad array marker.
    struct SimulationEnv se = {0};
    const char *names[] = {"A", "BB"};
    const size_t array_header_size =
        CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t);
    // Skip the scalar payload and the first three env arrays to land on the atom-name header.
    const size_t payload_offset = sizeof(CheckpointHeader) + sizeof(CheckpointEnvPayload) +
                                  (array_header_size + (size_t)2 * sizeof(double)) +
                                  (array_header_size + (size_t)1 * sizeof(double)) +
                                  (array_header_size + (size_t)2 * sizeof(bool));
    CheckpointStatus save_status;
    CheckpointStatus load_status;
    FILE *file;
    int ch;

    se.flavor = FLAVOR_MC;
    se.lattice_type = BCC;
    se.rand_seed = 123u;
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

    save_status = checkpoint_save(temp_checkpoint_path, NULL, &se, NULL);
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

    load_status = checkpoint_load(temp_checkpoint_path, NULL, &se, NULL);
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

void test_fill_atom_payload_copies_selected_atom_fields(void)
{
    // Copy a single atom into its compact on-disk payload representation.
    Atom atom = {0};
    CheckpointAtomPayload payload = {0};

    atom.type = 7u;
    atom.energy = -3.5;
    atom.lattice[0] = 1;
    atom.lattice[1] = 2;
    atom.lattice[2] = 3;
    atom.bsradius = 0.75;

    fill_atom_payload(&payload, &atom);

    // Mutate the source afterward so the test proves the payload kept the original value.
    atom.type = 9u;

    TEST_ASSERT_EQUAL_UINT_MESSAGE(7u, payload.type, "atom type should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(-3.5, payload.energy,
                                     "atom energy should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, payload.lattice_u, "lattice u should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, payload.lattice_v, "lattice v should copy into the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, payload.lattice_w, "lattice w should copy into the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(0.75, payload.bsradius,
                                     "bond-sphere radius should copy into the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(atom.type, payload.type,
                                       "payload should keep the copied atom type");
}

void test_apply_atom_payload_restores_selected_atom_fields(void)
{
    // Restore a single atom from its compact payload and verify the saved fields overwrite the old
    // ones.
    CheckpointAtomPayload payload = {
        .type = 4u,
        .energy = 1.25,
        .lattice_u = -1,
        .lattice_v = -2,
        .lattice_w = -3,
        .bsradius = 2.5,
    };
    Atom atom = {0};

    atom.type = 99u;
    atom.energy = -8.0;
    atom.lattice[0] = 10;
    atom.lattice[1] = 11;
    atom.lattice[2] = 12;
    atom.bsradius = 9.0;

    // Apply the payload to the live atom structure.
    apply_atom_payload(&payload, &atom);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(4u, atom.type, "atom type should restore from the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(99u, atom.type, "atom type should overwrite the old value");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(1.25, atom.energy,
                                     "atom energy should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(-1, atom.lattice[0], "lattice u should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(-2, atom.lattice[1], "lattice v should restore from the payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(-3, atom.lattice[2], "lattice w should restore from the payload");
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(2.5, atom.bsradius,
                                     "bond-sphere radius should restore from the payload");
}

void test_write_atom_array_serializes_atom_payloads(void)
{
    // Serialize a small atom array so the test can inspect the header and payload layout directly.
    Atom atom0 = {0};
    Atom atom1 = {0};
    Atom *atom_refs[] = {&atom0, &atom1};
    CheckpointAtomPayload expected_payloads[2] = {0};
    const size_t header_size = CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t);
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;

    atom0.type = 1u;
    atom0.energy = -1.5;
    atom0.lattice[0] = 0;
    atom0.lattice[1] = 1;
    atom0.lattice[2] = 2;
    atom0.bsradius = 0.5;

    atom1.type = 3u;
    atom1.energy = 4.25;
    atom1.lattice[0] = 4;
    atom1.lattice[1] = 5;
    atom1.lattice[2] = 6;
    atom1.bsradius = 1.25;

    fill_atom_payload(&expected_payloads[0], &atom0);
    fill_atom_payload(&expected_payloads[1], &atom1);

    write_atom_array((const Atom **)atom_refs, 2u, &payload, &payload_bytes);

    TEST_ASSERT_EQUAL_UINT_MESSAGE((unsigned int)(header_size + sizeof(expected_payloads)),
                                   payload_bytes,
                                   "atom array should report the full serialized size");
    // The first bytes are the shared array magic written before the typed header.
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&((const uint16_t){CHECKPOINT_ARRAY_MAGIC}), payload,
                                     sizeof(uint16_t),
                                     "atom array should begin with the array magic");
    TEST_ASSERT_EQUAL_UINT16_MESSAGE(CAF_ATOMS, *(uint16_t *)(payload + sizeof(uint16_t)),
                                     "atom array should encode the atom flag");
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(2u,
                                     *(uint32_t *)(payload + sizeof(uint16_t) + sizeof(uint16_t)),
                                     "atom array should encode the atom count");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(expected_payloads, payload + header_size,
                                     sizeof(expected_payloads),
                                     "atom payload bytes should match the compact atom layout");

    free(payload);
}

void test_read_atom_array_parses_atom_payloads(void)
{
    // Build a synthetic atom array payload and make sure the parser returns the right count and
    // bytes.
    enum { header_size = CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t) };
    CheckpointAtomPayload expected_payloads[2] = {
        {.type = 5u,
         .energy = -2.0,
         .lattice_u = 1,
         .lattice_v = 2,
         .lattice_w = 3,
         .bsradius = 0.25},
        {.type = 6u,
         .energy = 3.5,
         .lattice_u = 4,
         .lattice_v = 5,
         .lattice_w = 6,
         .bsradius = 0.75},
    };
    uint8_t payload[header_size + sizeof(expected_payloads)] = {0};
    size_t bytes_read = 0u;
    uint32_t out_n = 0u;
    CheckpointAtomPayload *out_arr = NULL;

    build_array_header(payload, CAF_ATOMS, 2u);
    memcpy(payload + header_size, expected_payloads, sizeof(expected_payloads));

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK,
                                  read_atom_array(payload, &bytes_read, &out_n, &out_arr),
                                  "valid atom array payload should parse successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(2u, out_n, "atom array should report the expected length");
    TEST_ASSERT_EQUAL_UINT_MESSAGE((unsigned int)(header_size + sizeof(expected_payloads)),
                                   (unsigned int)bytes_read,
                                   "atom array should report the consumed byte count");
    TEST_ASSERT_NOT_NULL_MESSAGE(out_arr, "atom array should allocate output storage");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(expected_payloads, out_arr, sizeof(expected_payloads),
                                     "atom array should round-trip the payload bytes");

    free(out_arr);
}

void test_read_atom_array_rejects_corrupted_magic(void)
{
    // Corrupt the array magic to confirm the parser stops before consuming bytes.
    uint8_t payload[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    size_t bytes_read = 123u;
    uint32_t out_n = 456u;
    CheckpointAtomPayload *out_arr = (CheckpointAtomPayload *)0x1;

    build_array_header(payload, CAF_ATOMS, 1u);
    payload[0] ^= 0xFFu;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR,
                                  read_atom_array(payload, &bytes_read, &out_n, &out_arr),
                                  "corrupted atom array magic should be rejected");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(123u, (unsigned int)bytes_read,
                                   "failed parse should not advance the byte count");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(456u, out_n, "failed parse should leave the length unchanged");
    TEST_ASSERT_EQUAL_PTR_MESSAGE((CheckpointAtomPayload *)0x1, out_arr,
                                  "failed parse should leave the output pointer untouched");
}

void test_validate_array_header_accepts_valid_header(void)
{
    // The raw header bytes should decode cleanly when the magic, flag, and length match.
    uint8_t header[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    uint16_t out_flag = 0u;
    uint32_t out_n = 0u;

    build_array_header(header, CAF_NN_ENERGY, 3u);

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, validate_array_header(header, &out_flag, &out_n),
                                  "valid header should parse successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(CAF_NN_ENERGY, out_flag,
                                   "valid header should report the expected flag");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(3u, out_n, "valid header should report the expected length");
}

void test_validate_array_header_rejects_invalid_magic(void)
{
    // Flip the magic bytes to make sure the header validator rejects a malformed record.
    uint8_t header[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    uint16_t out_flag = 0xFFFFu;
    uint32_t out_n = 0xFFFFFFFFu;

    build_array_header(header, CAF_IS_SOLUBLE, 2u);
    header[0] ^= 0xFFu;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR,
                                  validate_array_header(header, &out_flag, &out_n),
                                  "corrupted magic should be rejected");
}

void test_checkpoint_header_has_valid_checksum_accepts_matching_payload(void)
{
    // A checksum computed from the same payload should validate exactly.
    const uint8_t payload[] = {0x10u, 0x20u, 0x30u, 0x40u};
    CheckpointHeader header = {0};

    checkpoint_header_init(&header);
    checkpoint_header_finalize(&header, (uint32_t)sizeof(payload),
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

    checkpoint_header_init(&header);
    checkpoint_header_finalize(&header, (uint32_t)sizeof(payload),
                               checkpoint_checksum32(payload, sizeof(payload)));

    TEST_ASSERT_FALSE_MESSAGE(
        checkpoint_header_has_valid_checksum(&header, corrupted_payload, sizeof(corrupted_payload)),
        "corrupted payload should fail checksum validation");
}

void test_read_array_parses_double_array_payload(void)
{
    // Parse a numeric array and verify the helper reports the bytes it consumed.
    const double expected[] = {1.25, -2.5, 3.75};
    uint8_t payload[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t) +
                    sizeof(expected)] = {0};
    uint16_t out_flag = 0u;
    uint32_t out_n = 0u;
    size_t bytes_read = 0u;
    double *out_arr = NULL;

    build_array_header(payload, CAF_NN_ENERGY, (uint32_t)(sizeof(expected) / sizeof(expected[0])));
    memcpy(payload + CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t), expected,
           sizeof(expected));

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_OK, read_array(payload, &bytes_read, &out_flag, &out_n, (void **)&out_arr),
        "valid numeric array payload should parse successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(CAF_NN_ENERGY, out_flag,
                                   "parsed array should report the expected flag");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(3u, out_n, "parsed array should report the expected length");
    TEST_ASSERT_EQUAL_UINT_MESSAGE((unsigned int)(CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) +
                                                  sizeof(uint32_t) + sizeof(expected)),
                                   (unsigned int)bytes_read,
                                   "parsed array should report the consumed byte count");
    TEST_ASSERT_NOT_NULL_MESSAGE(out_arr, "parsed array should allocate output storage");
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY_MESSAGE(expected, out_arr, 3,
                                           "parsed array should round-trip element values");

    free(out_arr);
}

void test_read_array_rejects_corrupted_magic_without_consuming_bytes(void)
{
    // Corrupt the array magic and confirm the helper leaves the output parameters untouched.
    uint8_t payload[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    uint16_t out_flag = 0u;
    uint32_t out_n = 0u;
    size_t bytes_read = 123u;
    void *out_arr = (void *)0x1;

    build_array_header(payload, CAF_SUBSTRATE_COMPOSITION, 2u);
    payload[0] ^= 0xFFu;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR,
                                  read_array(payload, &bytes_read, &out_flag, &out_n, &out_arr),
                                  "corrupted array magic should be rejected");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(123u, (unsigned int)bytes_read,
                                   "failed parse should not advance the byte count");
    TEST_ASSERT_EQUAL_PTR_MESSAGE((void *)0x1, out_arr,
                                  "failed parse should leave the output pointer untouched");
}

void test_write_array_header_into_payload_writes_expected_header_bytes(void)
{
    // Write a standalone array header into a preallocated buffer and inspect the bytes directly.
    uint8_t buffer[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    uint8_t *buffer_ptr = buffer;
    const uint16_t expected_flag = CAF_ATOM_NAMES;
    const uint32_t expected_n = 5u;

    CheckpointStatus status =
        write_array_header_into_payload(expected_flag, expected_n, &buffer_ptr);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status,
                                  "valid output buffer should accept a header write");

    // Build an expected magic value in-place so the test can compare against raw bytes.
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&((const uint16_t){CHECKPOINT_ARRAY_MAGIC}), buffer,
                                     sizeof(uint16_t),
                                     "written header should start with the array magic");

    TEST_ASSERT_EQUAL_UINT16_MESSAGE(expected_flag, *(uint16_t *)(buffer + sizeof(uint16_t)),
                                     "written header should encode the array flag");
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(expected_n,
                                     *(uint32_t *)(buffer + sizeof(uint16_t) + sizeof(uint16_t)),
                                     "written header should encode the element count");
}

void test_write_array_header_into_payload_rejects_null_output_buffer(void)
{
    // A missing output buffer should fail immediately.
    uint8_t *buffer_ptr = NULL;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR,
                                  write_array_header_into_payload(CAF_ATOMS, 1u, &buffer_ptr),
                                  "null output storage should be rejected");
}

void test_write_array_serializes_double_array_payload(void)
{
    // Serialize a numeric array and compare the complete byte layout against the expected bytes.
    const double expected_values[] = {1.5, -2.0, 7.25};
    const size_t header_size = CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t);
    const size_t expected_size = header_size + sizeof(expected_values);
    uint8_t expected_payload[header_size + sizeof(expected_values)];
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;

    build_array_header(expected_payload, CAF_NN_ENERGY,
                       (uint32_t)(sizeof(expected_values) / sizeof(expected_values[0])));
    memcpy(expected_payload + header_size, expected_values, sizeof(expected_values));

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_OK,
        write_array(CAF_NN_ENERGY, (uint32_t)(sizeof(expected_values) / sizeof(expected_values[0])),
                    expected_values, sizeof(expected_values[0]), &payload, &payload_bytes),
        "valid array data should serialize successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE((unsigned int)expected_size, payload_bytes,
                                   "serialized array should report the full byte count");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(expected_payload, payload, expected_size,
                                     "serialized array should match the expected byte layout");

    free(payload);
}

void test_write_array_rejects_null_payload_destination(void)
{
    // The helper should reject a missing destination buffer before allocating anything.
    const double values[] = {3.0, 4.0};
    uint32_t payload_bytes = 0u;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_ERROR,
        write_array(CAF_NN_ENERGY, 2u, values, sizeof(values[0]), NULL, &payload_bytes),
        "null payload destination should be rejected");
}

void test_append_to_payload_appends_bytes_to_existing_payload(void)
{
    // Start with a tiny existing buffer, then append bytes and verify the old content stays intact.
    uint8_t *payload = (uint8_t *)malloc(2u);
    uint32_t payload_bytes = 2u;
    const uint8_t append_bytes[] = {0x30u, 0x40u, 0x50u};
    const uint8_t expected_bytes[] = {0x10u, 0x20u, 0x30u, 0x40u, 0x50u};

    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "payload should allocate");
    payload[0] = 0x10u;
    payload[1] = 0x20u;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_OK,
        append_to_payload(append_bytes, sizeof(append_bytes), &payload, &payload_bytes),
        "append should succeed with valid inputs");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(5u, payload_bytes,
                                   "payload byte count should include the appended bytes");
    TEST_ASSERT_EQUAL_UINT8_ARRAY_MESSAGE(expected_bytes, payload, 5u,
                                          "payload contents should preserve the original bytes");

    free(payload);
}

void test_append_to_payload_rejects_null_input_buffer(void)
{
    // Passing data without a source buffer should be treated as an immediate error.
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_ERROR, append_to_payload((const uint8_t[]){0xAAu}, 1u, &payload, &payload_bytes),
        "null payload input should be rejected");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(NULL, payload,
                                  "failed append should leave the payload pointer unchanged");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(0u, payload_bytes,
                                   "failed append should leave the payload size unchanged");
}

void test_write_array_magic_writes_expected_magic_bytes(void)
{
    // Verify the helper appends just the two-byte array magic and updates the byte count.
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;
    const uint16_t expected_magic = CHECKPOINT_ARRAY_MAGIC;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, write_array_magic(&payload, &payload_bytes),
                                  "writing array magic should succeed with valid inputs");
    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "array magic write should allocate payload bytes");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(2u, payload_bytes, "array magic write should append two bytes");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&expected_magic, payload, sizeof(expected_magic),
                                     "array magic bytes should match the expected value");

    free(payload);
}

void test_write_array_magic_rejects_null_length_pointer(void)
{
    // A null byte-count pointer should fail before any allocation or mutation occurs.
    uint8_t *payload = NULL;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, write_array_magic(&payload, NULL),
                                  "null byte-count storage should be rejected");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(NULL, payload,
                                  "failed write should leave the payload pointer unchanged");
}

void test_verify_payload_size_accepts_matching_file_size(void)
{
    // Build a synthetic file with the exact expected size so the verifier can seek to the payload.
    CheckpointHeader header = {0};
    FILE *file;
    const uint8_t payload[] = {0x11u, 0x22u, 0x33u, 0x44u};

    checkpoint_header_init(&header);
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

    checkpoint_header_init(&header);
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

void test_fill_env_array_payload_copies_selected_env_arrays(void)
{
    // Copy the env array pointers and string-length metadata into the temporary serialization
    // struct.
    struct SimulationEnv se = {0};
    CheckpointEnvArrPayload arr_payload = {0};
    double substrate_composition[] = {0.2, 0.8};
    double nn_energy[] = {1.25, -0.75, 0.5};
    bool is_soluble[] = {true, false};
    char name0[] = "Fe";
    char name1[] = "Ni";
    char *atom_names[] = {name0, name1};

    se.num_elements = 2;
    se.num_nn_types = 3;
    se.atom_names_cnt = 2;
    se.substrate_composition = substrate_composition;
    se.nn_energy = nn_energy;
    se.is_soluble = is_soluble;
    se.atom_names = atom_names;

    // Populate the payload view from the live env arrays.
    fill_env_array_payload(&arr_payload, &se);

    TEST_ASSERT_EQUAL_PTR_MESSAGE(substrate_composition, arr_payload.substrate_composition,
                                  "substrate composition pointer should be copied");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, arr_payload.n_substrate_composition,
                                  "substrate composition length should copy");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(nn_energy, arr_payload.nn_energy,
                                  "nn energy pointer should be copied");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, arr_payload.n_nn_energy, "nn energy length should copy");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(is_soluble, arr_payload.is_soluble,
                                  "solubility pointer should be copied");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, arr_payload.n_is_soluble, "solubility length should copy");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(atom_names, arr_payload.atom_names,
                                  "atom names pointer should be copied");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, arr_payload.n_atom_names, "atom name count should copy");
    TEST_ASSERT_NOT_NULL_MESSAGE(arr_payload.n_atom_names_str,
                                 "atom name length array should allocate");
    TEST_ASSERT_EQUAL_INT_MESSAGE((int)strlen(name0) + 1, arr_payload.n_atom_names_str[0],
                                  "first atom name length should include the terminator");
    TEST_ASSERT_EQUAL_INT_MESSAGE((int)strlen(name1) + 1, arr_payload.n_atom_names_str[1],
                                  "second atom name length should include the terminator");

    free(arr_payload.n_atom_names_str);
}

void test_apply_env_arrays_attaches_payload_pointers_to_env(void)
{
    // Attach the payload-owned arrays back to a blank env struct and verify the pointers move over.
    struct SimulationEnv se = {0};
    CheckpointEnvArrPayload arr_payload = {0};
    double substrate_composition[] = {0.3, 0.7};
    double nn_energy[] = {2.0};
    bool is_soluble[] = {false, true};
    char name0[] = "A";
    char name1[] = "B";
    char *atom_names[] = {name0, name1};

    arr_payload.substrate_composition = substrate_composition;
    arr_payload.nn_energy = nn_energy;
    arr_payload.is_soluble = is_soluble;
    arr_payload.atom_names = atom_names;
    arr_payload.n_atom_names = 2;

    // The helper should transfer ownership by copying the pointers into the env struct.
    apply_env_arrays(&arr_payload, &se);

    TEST_ASSERT_EQUAL_PTR_MESSAGE(substrate_composition, se.substrate_composition,
                                  "substrate composition pointer should attach to the env");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(nn_energy, se.nn_energy,
                                  "nn energy pointer should attach to the env");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(is_soluble, se.is_soluble,
                                  "solubility pointer should attach to the env");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(atom_names, se.atom_names,
                                  "atom name pointer should attach to the env");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, se.atom_names_cnt, "atom name count should attach to the env");
}

void test_write_and_read_env_arrays_round_trip(void)
{
    // Serialize the env arrays and parse them back to prove the standalone helpers agree on layout.
    struct SimulationEnv se = {0};
    CheckpointEnvArrPayload write_payload = {0};
    CheckpointEnvArrPayload read_payload = {0};
    double substrate_composition[] = {0.4, 0.6};
    double nn_energy[] = {1.0, 2.0};
    bool is_soluble[] = {true, false};
    char name0[] = "Al";
    char name1[] = "Cu";
    char *atom_names[] = {name0, name1};
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;
    size_t bytes_read = 0u;

    se.num_elements = 2;
    se.num_nn_types = 2;
    se.atom_names_cnt = 2;
    se.substrate_composition = substrate_composition;
    se.nn_energy = nn_energy;
    se.is_soluble = is_soluble;
    se.atom_names = atom_names;

    // Build the write-side view, serialize it, then free only the metadata allocated for names.
    fill_env_array_payload(&write_payload, &se);
    write_env_arrays(&write_payload, &payload, &payload_bytes);
    free(write_payload.n_atom_names_str);

    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "serialized env-array payload should allocate bytes");
    TEST_ASSERT_GREATER_THAN_UINT_MESSAGE(0u, payload_bytes,
                                          "serialized env-array payload should be non-empty");

    // Parse the raw bytes back into a fresh payload view.
    read_env_arrays(payload, &bytes_read, &read_payload);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(payload_bytes, (unsigned int)bytes_read,
                                   "read helper should consume the full payload");
    TEST_ASSERT_NOT_NULL_MESSAGE(read_payload.substrate_composition,
                                 "substrate composition should round-trip through serialization");
    TEST_ASSERT_NOT_NULL_MESSAGE(read_payload.nn_energy,
                                 "nn energy should round-trip through serialization");
    TEST_ASSERT_NOT_NULL_MESSAGE(read_payload.is_soluble,
                                 "solubility should round-trip through serialization");
    TEST_ASSERT_NOT_NULL_MESSAGE(read_payload.atom_names,
                                 "atom names should round-trip through serialization");
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY_MESSAGE(substrate_composition,
                                           read_payload.substrate_composition, 2,
                                           "substrate composition values should round-trip");
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY_MESSAGE(nn_energy, read_payload.nn_energy, 2,
                                           "nn energy values should round-trip");
    TEST_ASSERT_TRUE_MESSAGE(read_payload.is_soluble[0],
                             "first solubility value should round-trip");
    TEST_ASSERT_FALSE_MESSAGE(read_payload.is_soluble[1],
                              "second solubility value should round-trip");
    TEST_ASSERT_EQUAL_INT_MESSAGE(2, read_payload.n_atom_names,
                                  "atom name count should round-trip through serialization");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Al", read_payload.atom_names[0],
                                     "first atom name should round-trip");
    TEST_ASSERT_EQUAL_STRING_MESSAGE("Cu", read_payload.atom_names[1],
                                     "second atom name should round-trip");

    free(payload);
    free(read_payload.atom_names[0]);
    free(read_payload.atom_names[1]);
    free(read_payload.atom_names);
    free(read_payload.substrate_composition);
    free(read_payload.nn_energy);
    free(read_payload.is_soluble);
}

void test_read_env_arrays_rejects_corrupted_header(void)
{
    // Corrupt the serialized bytes before parsing so the helper should fail without touching
    // outputs.
    struct SimulationEnv se = {0};
    CheckpointEnvArrPayload write_payload = {0};
    CheckpointEnvArrPayload read_payload = {0};
    double substrate_composition[] = {0.4, 0.6};
    double nn_energy[] = {1.0};
    bool is_soluble[] = {true, false};
    char name0[] = "Al";
    char name1[] = "Cu";
    char *atom_names[] = {name0, name1};
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;
    size_t bytes_read = 123u;

    se.num_elements = 2;
    se.num_nn_types = 1;
    se.atom_names_cnt = 2;
    se.substrate_composition = substrate_composition;
    se.nn_energy = nn_energy;
    se.is_soluble = is_soluble;
    se.atom_names = atom_names;

    fill_env_array_payload(&write_payload, &se);
    write_env_arrays(&write_payload, &payload, &payload_bytes);
    free(write_payload.n_atom_names_str);

    // Flip the first byte of the array magic to simulate corruption.
    payload[0] ^= 0xFFu;

    read_env_arrays(payload, &bytes_read, &read_payload);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(123u, (unsigned int)bytes_read,
                                   "failed parse should leave the consumed-byte count unchanged");
    TEST_ASSERT_NULL_MESSAGE(read_payload.substrate_composition,
                             "failed parse should not attach substrate composition");
    TEST_ASSERT_NULL_MESSAGE(read_payload.atom_names, "failed parse should not attach atom names");

    free(payload);
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

    checkpoint_header_set_magic(&header);

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

    checkpoint_header_init(&header);

    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_is_valid(&header),
                             "initialized header should be valid");
}

void test_checkpoint_header_is_valid_rejects_invalid_structure(void)
{
    // Break the stored header size so the combined check rejects the structure.
    CheckpointHeader header = {0};

    checkpoint_header_init(&header);
    header.header_bytes = 0u;

    TEST_ASSERT_FALSE_MESSAGE(checkpoint_header_is_valid(&header),
                              "header with invalid structure should be rejected");
}

void test_checkpoint_header_has_valid_structure_accepts_initialized_header(void)
{
    // The initialized header should meet the version/size/reserved field requirements.
    CheckpointHeader header = {0};

    checkpoint_header_init(&header);

    TEST_ASSERT_TRUE_MESSAGE(checkpoint_header_has_valid_structure(&header),
                             "initialized header should have a valid structure");
}

void test_checkpoint_header_has_valid_structure_rejects_corrupted_header_bytes(void)
{
    // Zeroing the stored size should be enough to make the structure invalid.
    CheckpointHeader header = {0};

    checkpoint_header_init(&header);
    header.header_bytes = 0u;

    TEST_ASSERT_FALSE_MESSAGE(checkpoint_header_has_valid_structure(&header),
                              "corrupted header bytes should fail structure validation");
}

void test_checkpoint_header_init_and_finalize(void)
{
    // Initialize the header, then finalize it with an explicit payload size and checksum.
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
    // Write a header to disk and read it back unchanged.
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
    // Saving with no state or env should still produce a valid header-only file.
    CheckpointStatus save_status;
    CheckpointStatus load_status;

    save_status = checkpoint_save(temp_checkpoint_path, NULL, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "header-only checkpoint should save successfully");

    load_status = checkpoint_load(temp_checkpoint_path, NULL, NULL, NULL);
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

    save_status = checkpoint_save(temp_checkpoint_path, &ss, NULL, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "state scalar checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));

    load_status = checkpoint_load(temp_checkpoint_path, &ss, NULL, NULL);
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

    save_status = checkpoint_save(temp_checkpoint_path, NULL, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "env scalar checkpoint should save successfully");

    memset(&se, 0, sizeof(se));

    load_status = checkpoint_load(temp_checkpoint_path, NULL, &se, NULL);
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

    save_status = checkpoint_save(temp_checkpoint_path, &ss, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "combined state+env checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));
    memset(&se, 0, sizeof(se));

    load_status = checkpoint_load(temp_checkpoint_path, &ss, &se, NULL);
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
    // The restore path expects this env array to exist when rebuilding atom metadata.
    se.is_soluble = (bool *)malloc((size_t)se.num_elements * sizeof(bool));
    TEST_ASSERT_NOT_NULL_MESSAGE(se.is_soluble, "is_soluble array should allocate");
    for (int i = 0; i < se.num_elements; ++i) {
        se.is_soluble[i] = false;
    }

    // Point the state at the single atom so checkpoint_save writes the atom array.
    ss.atom_cnt = 1;
    ss.atom_arr = atom_refs;

    save_status = checkpoint_save(temp_checkpoint_path, &ss, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "atom checkpoint should save successfully");

    memset(&ss, 0, sizeof(ss));

    load_status = checkpoint_load(temp_checkpoint_path, &ss, &se, NULL);
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
    // Loading a missing path should fail cleanly rather than crashing or creating files.
    CheckpointStatus load_status;

    remove(temp_checkpoint_path);

    load_status = checkpoint_load(temp_checkpoint_path, NULL, NULL, NULL);
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

    save_status = checkpoint_save(temp_checkpoint_path, &ss, NULL, NULL);
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

    load_status = checkpoint_load(temp_checkpoint_path, &ss, NULL, NULL);
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
    save_status = checkpoint_save(temp_checkpoint_path, &ss_orig, &se, NULL);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, save_status,
                                  "checkpoint with atoms should save successfully");

    // Allocate the zone grid that checkpoint_load will fill while rebuilding derived data.
    initialize_zones(&ss_restored.zone_arr, &se);
    TEST_ASSERT_NOT_NULL(ss_restored.zone_arr);

    // Load the checkpoint and let the restore path rebuild the derived structures.
    load_status = checkpoint_load(temp_checkpoint_path, &ss_restored, &se, NULL);
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
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_checkpoint_checksum32_is_deterministic);
    RUN_TEST(test_fill_state_payload_copies_selected_state_scalars);
    RUN_TEST(test_apply_state_payload_to_simstate_restores_selected_state_scalars);
    RUN_TEST(test_fill_env_payload_copies_selected_env_scalars);
    RUN_TEST(test_apply_env_payload_to_config_restores_selected_env_scalars);
    RUN_TEST(test_checkpoint_save_and_load_env_arrays_and_atom_names);
    RUN_TEST(test_checkpoint_load_rejects_corrupted_env_atom_names_header);
    RUN_TEST(test_fill_atom_payload_copies_selected_atom_fields);
    RUN_TEST(test_apply_atom_payload_restores_selected_atom_fields);
    RUN_TEST(test_write_atom_array_serializes_atom_payloads);
    RUN_TEST(test_read_atom_array_parses_atom_payloads);
    RUN_TEST(test_read_atom_array_rejects_corrupted_magic);
    RUN_TEST(test_validate_array_header_accepts_valid_header);
    RUN_TEST(test_validate_array_header_rejects_invalid_magic);
    RUN_TEST(test_checkpoint_header_has_valid_checksum_accepts_matching_payload);
    RUN_TEST(test_checkpoint_header_has_valid_checksum_rejects_corrupted_payload);
    RUN_TEST(test_read_array_parses_double_array_payload);
    RUN_TEST(test_read_array_rejects_corrupted_magic_without_consuming_bytes);
    RUN_TEST(test_write_array_header_into_payload_writes_expected_header_bytes);
    RUN_TEST(test_write_array_header_into_payload_rejects_null_output_buffer);
    RUN_TEST(test_write_array_serializes_double_array_payload);
    RUN_TEST(test_write_array_rejects_null_payload_destination);
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
