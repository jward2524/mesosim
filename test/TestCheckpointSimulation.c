#include "CheckpointSimulation.h"
#include "TUtils.h"
#include "unity.h"
#include <stdlib.h>
#include <string.h>

// TODO: set up ss and se in setUp so it cleans on test fail
void setUp(void)
{
}

void tearDown(void)
{
}

void test_fill_state_payload_copies_selected_state_scalars(void)
{
    // Seed the state with values that should survive the copy into the compact payload.
    struct SimulationState ss = {0};
    SimStatePayload payload = {0};

    ss.iter = 17ul;
    ss.mmc_step = 23ul;
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

    pack_simstate(&payload, &ss);

    // Mutate the source after filling the payload to confirm the payload owns its copied values.
    ss.iter = 999ul;

    TEST_ASSERT_EQUAL_UINT_MESSAGE(17u, payload.iter, "iter should copy into the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(ss.iter, payload.iter,
                                       "payload should keep the copied iter value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(23u, payload.mmc_step, "mmc steps should copy into the payload");
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
    SimStatePayload payload = {
        .iter = 11ul,
        .mmc_step = 22ul,
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
        .mmc_step = 666ul,
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
    unpack_state_payload_to_simstate(&payload, &ss);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(11u, ss.iter, "iter should restore from the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(555u, ss.iter, "iter should overwrite the old value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(22u, ss.mmc_step, "mmc steps should restore from payload");
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
    SimEnvPayload payload = {0};

    se.flavor = FLAVOR_KMC;
    se.rand_state = (RandomState){.u = 12345, .v = 67890, .w = 24680};
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

    pack_simenv(&payload, &se);

    // Change the source after copying so the payload value is the one under test.
    se.flavor = FLAVOR_MC;

    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, payload.flavor,
                                   "flavor should copy into the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(se.flavor, payload.flavor,
                                       "payload should keep the copied flavor value");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(12345, payload.rand_state.u,
                                   "random state u should copy into the payload");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(67890, payload.rand_state.v,
                                   "random state v should copy into the payload");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(24680, payload.rand_state.w,
                                   "random state w should copy into the payload");
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
    SimEnvPayload payload = {
        .flavor = FLAVOR_MC,
        .rand_state = (RandomState){.u = 12345, .v = 67890, .w = 24680},
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
    // config should reflect the values in payload if loaded correctly
    struct SimulationConfig config = {
        .flavor = FLAVOR_KMC,
        .rand_state = (RandomState){.u = 12345, .v = 67890, .w = 24680},
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
    unpack_simenv_payload_to_config(&payload, &config);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(FLAVOR_MC, config.flavor,
                                   "flavor should restore from the payload");
    TEST_ASSERT_NOT_EQUAL_UINT_MESSAGE(FLAVOR_KMC, config.flavor,
                                       "flavor should overwrite the old value");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&payload.rand_state, &config.rand_state, sizeof(RandomState),
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
void test_fill_atom_payload_copies_selected_atom_fields(void)
{
    // Copy a single atom into its compact on-disk payload representation.
    Atom atom = {0};
    AtomPayload payload = {0};

    atom.type = 7u;
    atom.energy = -3.5;
    atom.lattice[0] = 1;
    atom.lattice[1] = 2;
    atom.lattice[2] = 3;
    atom.bsradius = 0.75;

    pack_atom(&payload, &atom);

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
    AtomPayload payload = {
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
    unpack_atom(&payload, &atom);

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
    AtomPayload expected_payloads[2] = {0};
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

    pack_atom(&expected_payloads[0], &atom0);
    pack_atom(&expected_payloads[1], &atom1);

    serialize_atom_array((const Atom **)atom_refs, 2u, &payload, &payload_bytes);

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
    AtomPayload expected_payloads[2] = {
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
    AtomPayload *out_arr = NULL;

    build_array_header(payload, CAF_ATOMS, 2u);
    memcpy(payload + header_size, expected_payloads, sizeof(expected_payloads));

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK,
                                  unserialize_atom_array(payload, &bytes_read, &out_n, &out_arr),
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
    AtomPayload *out_arr = (AtomPayload *)0x1;

    build_array_header(payload, CAF_ATOMS, 1u);
    payload[0] ^= 0xFFu;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR,
                                  unserialize_atom_array(payload, &bytes_read, &out_n, &out_arr),
                                  "corrupted atom array magic should be rejected");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(123u, (unsigned int)bytes_read,
                                   "failed parse should not advance the byte count");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(456u, out_n, "failed parse should leave the length unchanged");
    TEST_ASSERT_EQUAL_PTR_MESSAGE((AtomPayload *)0x1, out_arr,
                                  "failed parse should leave the output pointer untouched");
}
void test_fill_env_array_payload_copies_selected_env_arrays(void)
{
    // Copy the env array pointers and string-length metadata into the temporary serialization
    // struct.
    struct SimulationEnv se = {0};
    SimEnvArrPayload arr_payload = {0};
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
    pack_env_arrays(&arr_payload, &se);

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
    SimEnvArrPayload arr_payload = {0};
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
    unpack_env_arrays(&arr_payload, &se);

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
    SimEnvArrPayload write_payload = {0};
    SimEnvArrPayload read_payload = {0};
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
    pack_env_arrays(&write_payload, &se);
    serialize_env_arrays(&write_payload, &payload, &payload_bytes);
    free(write_payload.n_atom_names_str);

    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "serialized env-array payload should allocate bytes");
    TEST_ASSERT_GREATER_THAN_UINT_MESSAGE(0u, payload_bytes,
                                          "serialized env-array payload should be non-empty");

    // Parse the raw bytes back into a fresh payload view.
    unserialize_env_arrays(payload, &bytes_read, &read_payload);

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
    SimEnvArrPayload write_payload = {0};
    SimEnvArrPayload read_payload = {0};
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

    pack_env_arrays(&write_payload, &se);
    serialize_env_arrays(&write_payload, &payload, &payload_bytes);
    free(write_payload.n_atom_names_str);

    // Flip the first byte of the array magic to simulate corruption.
    payload[0] ^= 0xFFu;

    unserialize_env_arrays(payload, &bytes_read, &read_payload);

    TEST_ASSERT_EQUAL_UINT_MESSAGE(123u, (unsigned int)bytes_read,
                                   "failed parse should leave the consumed-byte count unchanged");
    TEST_ASSERT_NULL_MESSAGE(read_payload.substrate_composition,
                             "failed parse should not attach substrate composition");
    TEST_ASSERT_NULL_MESSAGE(read_payload.atom_names, "failed parse should not attach atom names");

    free(payload);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_fill_state_payload_copies_selected_state_scalars);
    RUN_TEST(test_apply_state_payload_to_simstate_restores_selected_state_scalars);

    RUN_TEST(test_fill_env_payload_copies_selected_env_scalars);
    RUN_TEST(test_apply_env_payload_to_config_restores_selected_env_scalars);

    RUN_TEST(test_fill_atom_payload_copies_selected_atom_fields);
    RUN_TEST(test_apply_atom_payload_restores_selected_atom_fields);
    RUN_TEST(test_write_atom_array_serializes_atom_payloads);
    RUN_TEST(test_read_atom_array_parses_atom_payloads);
    RUN_TEST(test_read_atom_array_rejects_corrupted_magic);

    UNITY_END();
    return 0;
}
