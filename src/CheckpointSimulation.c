#include "CheckpointSimulation.h"
#include "Atoms.h"
#include "State.h"
#include "Utils.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>

// copy the SimulationState fields we care about into a compact on-disk payload.
void pack_simstate(SimStatePayload *payload, const struct SimulationState *ss)
{
    payload->iter = ss->iter;
    payload->mmc_steps = ss->mmc_steps;
    payload->final_iteration = ss->final_iteration;
    payload->run_stime = ss->run_stime;
    payload->simulation_should_kill_itself = ss->simulation_should_kill_itself;
    payload->elapsed_stime = ss->elapsed_stime;
    payload->sim_end_type = ss->sim_end_type;
    payload->frequency_sum = ss->frequency_sum;
    payload->total_internal_energy = ss->total_internal_energy;
    payload->temperature = ss->temperature;
    payload->overpotential = ss->overpotential;
    payload->total_atoms_dissolved = ss->total_atoms_dissolved;
}

// copy the serialized payload back into a SimulationState so a restored run can resume with the
// same counters and scalar conditions.
void unpack_state_payload_to_simstate(const SimStatePayload *payload, struct SimulationState *ss)
{
    ss->iter = payload->iter;
    ss->mmc_steps = payload->mmc_steps;
    ss->final_iteration = payload->final_iteration;
    ss->run_stime = payload->run_stime;
    ss->simulation_should_kill_itself = payload->simulation_should_kill_itself; // XXX
    ss->elapsed_stime = payload->elapsed_stime;
    ss->sim_end_type = payload->sim_end_type;
    ss->frequency_sum = payload->frequency_sum; // XXX
    ss->temperature = payload->temperature;
    ss->overpotential = payload->overpotential;
    ss->total_atoms_dissolved = payload->total_atoms_dissolved;

    // XXX: updated by refresh_transitions
    ss->total_internal_energy = payload->total_internal_energy;
}

// copy the immutable SimulationEnv fields into a compact on-disk payload for restore validation.
void pack_simenv(SimEnvPayload *payload, const struct SimulationEnv *se)
{
    payload->flavor = se->flavor;
    payload->rand_seed = se->rand_seed;

    payload->system_size_x = se->system_size_x;
    payload->system_size_y = se->system_size_y;
    payload->system_size_z = se->system_size_z;

    payload->num_elements = se->num_elements;
    // TODO: change these derived fields to be re-calculated instead of passed through
    // requires functionalization/reuse of code in finalize_nne
    payload->num_nn_levels = se->num_nn_levels;
    payload->num_bond_types = se->num_bond_types;
    payload->num_nn_types = se->num_nn_types;
    // payload->num_neighbor_types = se->num_neighbor_types;
    payload->num_transition_vectors = se->num_transition_vectors;

    payload->overpotential_ramp_rate = se->overpotential_ramp_rate;
    payload->max_overpotential = se->max_overpotential;

    payload->dissolution = se->dissolution;
    payload->atom_names_cnt = se->atom_names_cnt;

    // TODO: lattice type is only in se for this transfer
    // is there a better way?
    payload->lattice_type = se->lattice_type;
}

// restore the immutable SimulationEnv config from a saved payload.
void unpack_simenv_payload_to_config(const SimEnvPayload *payload, struct SimulationConfig *config)
{
    config->flavor = payload->flavor;
    config->rand_seed = payload->rand_seed;
    config->system_size_x = payload->system_size_x;
    config->system_size_y = payload->system_size_y;
    config->system_size_z = payload->system_size_z;
    config->num_elements = payload->num_elements;
    config->num_nn_levels = payload->num_nn_levels;
    config->num_bond_types = payload->num_bond_types;

    // TODO: this might be handled elsewhere
    config->num_nn_types = payload->num_nn_types;

    // se->num_neighbor_types = payload->num_neighbor_types;
    // config->num_transition_vectors = payload->num_transition_vectors;
    config->overpotential_ramp_rate = payload->overpotential_ramp_rate;
    config->max_overpotential = payload->max_overpotential;
    config->dissolution = payload->dissolution;
    config->atom_names_cnt = payload->atom_names_cnt;
    config->lattice_type = payload->lattice_type;
}

void pack_env_arrays(SimEnvArrPayload *arr_payload, const struct SimulationEnv *se)
{
    arr_payload->substrate_composition = se->substrate_composition;
    arr_payload->n_substrate_composition = (int)se->num_elements;

    arr_payload->nn_energy = se->nn_energy;
    arr_payload->n_nn_energy = (int)se->num_nn_types;

    arr_payload->is_soluble = se->is_soluble;
    arr_payload->n_is_soluble = (int)se->num_elements;

    arr_payload->atom_names = se->atom_names;
    arr_payload->n_atom_names = (int)se->atom_names_cnt;
    arr_payload->n_atom_names_str = (int *)malloc((size_t)arr_payload->n_atom_names * sizeof(int));
    for (int i = 0; i < arr_payload->n_atom_names; ++i) {
        arr_payload->n_atom_names_str[i] =
            (int)strlen(se->atom_names[i]) + 1; // include null terminator
    }
}

void serialize_env_arrays(const SimEnvArrPayload *arr_payload, uint8_t **p_payload,
                          uint32_t *p_payload_bytes)
{
    CheckpointStatus status;
    status = serialize_array(
        (uint16_t)CAF_SUBSTRATE_COMPOSITION, (uint32_t)arr_payload->n_substrate_composition,
        arr_payload->substrate_composition, sizeof(double), p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write substrate composition array\n");
        return;
    }

    status = serialize_array((uint16_t)CAF_NN_ENERGY, (uint32_t)arr_payload->n_nn_energy,
                             arr_payload->nn_energy, sizeof(double), p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write nn energy array\n");
        return;
    }
    status = serialize_array((uint16_t)CAF_IS_SOLUBLE, (uint32_t)arr_payload->n_is_soluble,
                             arr_payload->is_soluble, sizeof(bool), p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write is_soluble array\n");
        return;
    }

    // write atom names array header first
    // use arr=NULL and elem_size=0 to write the strings separately
    status = serialize_array((uint16_t)CAF_ATOM_NAMES, (uint32_t)arr_payload->n_atom_names, NULL, 0,
                             p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write atom names array header\n");
    }

    // for atom names, we need to write each string separately since they may have different
    // lengths
    for (int i = 0; i < arr_payload->n_atom_names; i++) {
        const char *name = arr_payload->atom_names[i];
        uint32_t name_len = (uint32_t)strlen(name) + 1; // include null terminator
        status = serialize_array(CAF_ATOM_NAME_STR, name_len, name, sizeof(char), p_payload,
                                 p_payload_bytes);
        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to write atom name string for %s\n", name);
            return;
        }
    }
}

void unserialize_env_arrays(uint8_t *payload, size_t *total_bytes_read,
                            SimEnvArrPayload *arr_payload)
{
    CheckpointStatus status;
    uint16_t flag = 0;
    uint32_t n = 0;
    void *arr = NULL;
    size_t bytes_read = 0;
    uint8_t *payload_ptr = payload;

    // the order of these reads must match the order of writes in fill_env_arrays

    // substrate composition
    status = unserialize_array(payload_ptr, &bytes_read, &flag, &n, &arr);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to read substrate composition array\n");
        return;
    }
    if (flag != CAF_SUBSTRATE_COMPOSITION) {
        fprintf(stderr, "Checkpoint error: expected substrate composition array flag %u, got %u\n",
                CAF_SUBSTRATE_COMPOSITION, flag);
        free(arr);
        return;
    }
    payload_ptr += bytes_read;
    arr_payload->substrate_composition = (double *)arr;

    // nn energy
    status = unserialize_array(payload_ptr, &bytes_read, &flag, &n, &arr);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to read nn energy array\n");
        return;
    }
    if (flag != CAF_NN_ENERGY) {
        fprintf(stderr, "Checkpoint error: expected nn energy array flag %u, got %u\n",
                CAF_NN_ENERGY, flag);
        free(arr);
        return;
    }
    payload_ptr += bytes_read;
    arr_payload->nn_energy = (double *)arr;

    // is_soluble
    status = unserialize_array(payload_ptr, &bytes_read, &flag, &n, &arr);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to read is_soluble array\n");
        return;
    }
    if (flag != CAF_IS_SOLUBLE) {
        fprintf(stderr, "Checkpoint error: expected is_soluble array flag %u, got %u\n",
                CAF_IS_SOLUBLE, flag);
        free(arr);
        return;
    }
    payload_ptr += bytes_read;
    arr_payload->is_soluble = (bool *)arr;

    // atom names array header
    status = unserialize_array(payload_ptr, &bytes_read, &flag, &n, &arr);
    // arr is still unset if flag=CAF_ATOM_NAMES
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to read atom names array\n");
        return;
    }
    if (flag != CAF_ATOM_NAMES) {
        fprintf(stderr, "Checkpoint error: expected atom names array flag %u, got %u\n",
                CAF_ATOM_NAMES, flag);
        return;
    }
    payload_ptr += bytes_read;
    arr_payload->n_atom_names = (int)n;
    arr_payload->atom_names = (char **)malloc(n * sizeof(char *));
    if (arr_payload->atom_names == NULL) {
        fprintf(stderr, "Checkpoint error: failed to allocate memory for atom names: %s\n",
                strerror(errno));
        free(arr);
        return;
    }
    for (uint32_t i = 0; i < n; i++) {
        uint16_t name_flag;
        uint32_t name_n;
        char *name_arr;
        status =
            unserialize_array(payload_ptr, &bytes_read, &name_flag, &name_n, (void **)&name_arr);
        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to read atom name string array\n");
            free(arr);
            free(arr_payload->atom_names);
            return;
        }
        if (name_flag != CAF_ATOM_NAME_STR) {
            fprintf(stderr, "Checkpoint error: expected atom name string flag %u, got %u\n",
                    CAF_ATOM_NAME_STR, name_flag);
            free(arr);
            free(arr_payload->atom_names);
            return;
        }
        payload_ptr += bytes_read;
        arr_payload->atom_names[i] = name_arr;
    }
    *total_bytes_read = (size_t)(payload_ptr - payload);
}

void unpack_env_arrays(SimEnvArrPayload *arr_payload, struct SimulationEnv *se)
{
    // copying the malloc'd pointers in arr_payload to SimulationEnv

    se->substrate_composition = arr_payload->substrate_composition;
    // se->atoms_per_nn_level = arr_payload->atoms_per_nn_level;
    se->nn_energy = arr_payload->nn_energy;
    se->is_soluble = arr_payload->is_soluble;
    se->atom_names = arr_payload->atom_names;

    // se->num_nn_types = arr_payload->n_nn_energy;
    se->atom_names_cnt = arr_payload->n_atom_names;
}

// serialize a single atom into its checkpoint representation.
void pack_atom(AtomPayload *payload, const Atom *atom)
{
    payload->type = atom->type;
    payload->energy = atom->energy;
    payload->lattice_u = atom->lattice[0];
    payload->lattice_v = atom->lattice[1];
    payload->lattice_w = atom->lattice[2];
    payload->bsradius = atom->bsradius;
}

// restore a single atom from its checkpoint representation.
// note: transition_indices, neighbor_atom_idxs, and linked-list fields are zeroed;
// they will be rebuilt from the simulation state during restore.
void unpack_atom(const AtomPayload *payload, Atom *atom)
{
    atom->type = payload->type;
    atom->energy = payload->energy;
    atom->lattice[0] = payload->lattice_u;
    atom->lattice[1] = payload->lattice_v;
    atom->lattice[2] = payload->lattice_w;
    atom->bsradius = payload->bsradius;
}

void serialize_atom_array(const Atom **atom_arr, uint32_t n_atoms, uint8_t **p_payload,
                          uint32_t *p_payload_bytes)
{
    // write the atom array header first
    serialize_array((uint16_t)CAF_ATOMS, n_atoms, NULL, 0, p_payload, p_payload_bytes);

    // write each atom's data as a compact payload
    for (uint32_t i = 0; i < n_atoms; ++i) {
        AtomPayload atom_payload;
        pack_atom(&atom_payload, atom_arr[i]);
        append_to_payload(&atom_payload, sizeof(AtomPayload), p_payload, p_payload_bytes);
    }
}

CheckpointStatus unserialize_atom_array(const uint8_t *payload, size_t *total_bytes_read,
                                        uint32_t *out_n, AtomPayload **out_atom_payload_arr)
{
    CheckpointStatus status;
    uint16_t flag = 0;
    uint32_t n = 0;

    status = unserialize_array(payload, total_bytes_read, &flag, &n, (void **)out_atom_payload_arr);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to read atom array header\n");
        return status;
    }
    if (flag != CAF_ATOMS) {
        fprintf(stderr, "Checkpoint error: expected atom array flag %u, got %u\n", CAF_ATOMS, flag);
        free(out_atom_payload_arr);
        return CHECKPOINT_ERROR;
    }
    *out_n = n;
    return CHECKPOINT_OK;
}

CheckpointStatus unpack_atom_array(const AtomPayload *atom_payload_arr, uint32_t n_atoms,
                                   int num_elements, struct SimulationState *ss,
                                   struct SimulationEnv *se)
{
    for (uint32_t i = 0; i < n_atoms; ++i) {
        // ENHANCE: replace with add_atom?
        Atom *atom = (Atom *)malloc(sizeof(Atom));
        if (atom == NULL) {
            fprintf(stderr, "Checkpoint error: failed to allocate memory for atom %d: %s\n", i,
                    strerror(errno));
            // cleanup previously allocated atoms
            for (uint32_t j = 0; j < i; ++j) {
                free(ss->atom_arr[j]);
            }
            return CHECKPOINT_ERROR;
        }

        unpack_atom(&atom_payload_arr[i], atom);
        if (atom->type >= num_elements) {
            fprintf(stderr, "Checkpoint error: invalid atom type %u for atom %d\n", atom->type, i);
            free(atom);
            // cleanup previously allocated atoms
            for (uint32_t j = 0; j < i; ++j) {
                free(ss->atom_arr[j]);
            }
            return CHECKPOINT_ERROR;
        }

        // needs to be added to atom array early because other fxns rely on atom index
        ss->atom_arr[i] = atom;

        lattice2cartesian(atom->lattice, se->primitive_basis, atom->cartesian);

        // find which zone this atom belongs to
        int zone_u, zone_v, zone_w;
        findzone(&zone_u, &zone_v, &zone_w, atom->lattice[0], atom->lattice[1], atom->lattice[2],
                 se);
        add_atom_to_zone_list((long)i, zone_u, zone_v, zone_w, ss);

        // initialize members not stored in payload to default values
        atom->next_atom = -1;
        atom->previous_atom = -1;
        for (int j = 0; j < se->num_transition_vectors; ++j) {
            atom->neighbor_atom_idxs[j] = -1;
            atom->transition_indices[j] = -1;
        }
        atom->transition_indices[se->num_transition_vectors + se->dissolution] = -1;
    }
    ss->atom_cnt = (long)n_atoms;
    return CHECKPOINT_OK;
}
