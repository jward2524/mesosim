#ifndef CHECKPOINT_SIMULATION_H
#define CHECKPOINT_SIMULATION_H

#include "CheckpointUtils.h"
#include "State.h"
#include <stdbool.h>
#include <stdint.h>

#pragma pack(push, 1)
typedef struct {
    unsigned long iter;
    unsigned long mmc_steps;
    unsigned long final_iteration;
    double run_stime;
    bool simulation_should_kill_itself;
    double elapsed_stime;
    int sim_end_type;
    double frequency_sum;
    double total_internal_energy;
    double temperature;
    double overpotential;
    int total_atoms_dissolved;
} SimStatePayload;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct {
    unsigned flavor;
    unsigned rand_seed;

    double overpotential_ramp_rate;
    double max_overpotential;

    int system_size_x;
    int system_size_y;
    int system_size_z;

    int num_elements;
    int num_nn_levels;
    int num_bond_types;
    int num_nn_types;
    int num_transition_vectors;

    int dissolution;
    int atom_names_cnt;
    int lattice_type;
} SimEnvPayload;
#pragma pack(pop)

typedef struct SEAP {
    double *substrate_composition;
    int n_substrate_composition;
    double *nn_energy;
    int n_nn_energy;
    bool *is_soluble;
    int n_is_soluble;
    char **atom_names;
    int n_atom_names;
    int *n_atom_names_str; // size of each string in atom_names
} SimEnvArrPayload;

#pragma pack(push, 1)
typedef struct {
    unsigned char type;
    double energy;
    int lattice_u;
    int lattice_v;
    int lattice_w;
    double bsradius;
} AtomPayload;
#pragma pack(pop)

void pack_simstate(SimStatePayload *payload, const struct SimulationState *ss);
void unpack_state_payload_to_simstate(const SimStatePayload *payload, struct SimulationState *ss);

void pack_simenv(SimEnvPayload *payload, const struct SimulationEnv *se);
void unpack_simenv_payload_to_config(const SimEnvPayload *payload, struct SimulationConfig *config);
void pack_env_arrays(SimEnvArrPayload *arr_payload, const struct SimulationEnv *se);
void unpack_env_arrays(SimEnvArrPayload *arr_payload, struct SimulationEnv *se);
void serialize_env_arrays(const SimEnvArrPayload *arr_payload, uint8_t **p_payload,
                          uint32_t *p_payload_bytes);
void unserialize_env_arrays(uint8_t *payload, size_t *total_bytes_read,
                            SimEnvArrPayload *arr_payload);

void pack_atom(AtomPayload *payload, const Atom *atom);
void unpack_atom(const AtomPayload *payload, Atom *atom);
void serialize_atom_array(const Atom **atom_arr, uint32_t n_atoms, uint8_t **p_payload,
                          uint32_t *p_payload_bytes);
CheckpointStatus unserialize_atom_array(const uint8_t *payload, size_t *total_bytes_read,
                                        uint32_t *out_n, AtomPayload **out_atom_payload_arr);

// packing atom array is done in serialize_atom_array
CheckpointStatus unpack_atom_array(const AtomPayload *atom_payload_arr, uint32_t n_atoms,
                                   int num_elements, struct SimulationState *ss,
                                   struct SimulationEnv *se);

#endif // CHECKPOINT_SIMULATION_H
