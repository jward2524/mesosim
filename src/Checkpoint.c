#include "Checkpoint.h"
#include "Atoms.h"
#include "Initialization.h"
#include "Simulation.h"
#include "Utils.h"
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// copy the SimulationState fields we care about into a compact on-disk payload.
void fill_state_payload(CheckpointStatePayload *payload, const struct SimulationState *ss)
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
void apply_state_payload_to_simstate(const CheckpointStatePayload *payload,
                                     struct SimulationState *ss)
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
void fill_env_payload(CheckpointEnvPayload *payload, const struct SimulationEnv *se)
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
void apply_env_payload_to_config(const CheckpointEnvPayload *payload,
                                 struct SimulationConfig *config)
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

void fill_env_array_payload(CheckpointEnvArrPayload *arr_payload, const struct SimulationEnv *se)
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

CheckpointStatus append_to_payload(const void *data, size_t data_size, uint8_t **p_payload,
                                   uint32_t *p_payload_bytes)
{
    if (p_payload == NULL || p_payload_bytes == NULL || (data == NULL && data_size > 0)) {
        fprintf(stderr, "Checkpoint error: invalid output buffer for payload append\n");
        return CHECKPOINT_ERROR;
    }
    *p_payload = realloc(*p_payload, (size_t)*p_payload_bytes + data_size);
    if (*p_payload == NULL) {
        fprintf(stderr, "Checkpoint error: failed to allocate memory for payload: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }
    uint8_t *payload_end = *p_payload + *p_payload_bytes;
    memcpy(payload_end, data, data_size);
    *p_payload_bytes += (uint32_t)data_size;
    return CHECKPOINT_OK;
}

CheckpointStatus write_array_magic(uint8_t **p_payload, uint32_t *p_payload_bytes)
{
    const uint16_t arr_magic = CHECKPOINT_ARRAY_MAGIC;
    CheckpointStatus status =
        append_to_payload(&arr_magic, sizeof(arr_magic), p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write array magic: %s\n", strerror(errno));
    }
    return status;
}

// TODO: make type of flag smaller and make sure type of n is large enough for everything
CheckpointStatus write_array_header_into_payload(uint16_t flag, uint32_t n, uint8_t **p_payload)
{
    if (p_payload == NULL || *p_payload == NULL) {
        fprintf(stderr, "Checkpoint error: invalid output buffer for array header\n");
        return CHECKPOINT_ERROR;
    }

    const uint16_t arr_magic = CHECKPOINT_ARRAY_MAGIC;
    size_t header_size = sizeof(arr_magic) + sizeof(flag) + sizeof(n);
    uint8_t *header = malloc(header_size);
    if (header == NULL) {
        fprintf(stderr, "Checkpoint error: failed to allocate memory for array header: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    memcpy(header, &arr_magic, sizeof(arr_magic));
    memcpy(header + sizeof(arr_magic), &flag, sizeof(flag));
    memcpy(header + sizeof(arr_magic) + sizeof(flag), &n, sizeof(n));

    memcpy(*p_payload, header, header_size);

    free(header);
    return CHECKPOINT_OK;
}

/**
 * @brief Writes an array to the checkpoint payload with a header containing a magic number, flag,
 * and length. Will write only the array header if arr is NULL.
 *
 * @param flag
 * @param n
 * @param arr
 * @param elem_size
 * @param payload
 * @param payload_bytes
 * @return CheckpointStatus
 */
CheckpointStatus write_array(uint16_t flag, uint32_t n, const void *arr, size_t elem_size,
                             uint8_t **p_payload, uint32_t *p_payload_bytes)
{
    CheckpointStatus status;
    if (p_payload == NULL || p_payload_bytes == NULL) {
        fprintf(stderr, "Checkpoint error: invalid output buffer for array payload\n");
        return CHECKPOINT_ERROR;
    }

    const uint16_t arr_magic = CHECKPOINT_ARRAY_MAGIC;
    size_t header_size = sizeof(arr_magic) + sizeof(flag) + sizeof(n);
    size_t intermediate_size = header_size + elem_size * n;
    uint8_t *intermediate_payload = malloc(intermediate_size);
    if (intermediate_payload == NULL) {
        fprintf(stderr, "Checkpoint error: failed to allocate memory for array payload: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    status = write_array_header_into_payload(flag, n, &intermediate_payload);
    if (status != CHECKPOINT_OK) {
        free(intermediate_payload);
        return status;
    }
    if (arr != NULL) {
        memcpy(intermediate_payload + header_size, arr, elem_size * n);
    }

    status = append_to_payload(intermediate_payload, intermediate_size, p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write array to payload: %s\n",
                strerror(errno));
    }
    free(intermediate_payload);
    return status;
}

void write_env_arrays(const CheckpointEnvArrPayload *arr_payload, uint8_t **p_payload,
                      uint32_t *p_payload_bytes)
{
    CheckpointStatus status;
    status = write_array(
        (uint16_t)CAF_SUBSTRATE_COMPOSITION, (uint32_t)arr_payload->n_substrate_composition,
        arr_payload->substrate_composition, sizeof(double), p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write substrate composition array\n");
        return;
    }

    status = write_array((uint16_t)CAF_NN_ENERGY, (uint32_t)arr_payload->n_nn_energy,
                         arr_payload->nn_energy, sizeof(double), p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write nn energy array\n");
        return;
    }
    status = write_array((uint16_t)CAF_IS_SOLUBLE, (uint32_t)arr_payload->n_is_soluble,
                         arr_payload->is_soluble, sizeof(bool), p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write is_soluble array\n");
        return;
    }

    // write atom names array header first
    // use arr=NULL and elem_size=0 to write the strings separately
    status = write_array((uint16_t)CAF_ATOM_NAMES, (uint32_t)arr_payload->n_atom_names, NULL, 0,
                         p_payload, p_payload_bytes);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write atom names array header\n");
    }

    // for atom names, we need to write each string separately since they may have different
    // lengths
    for (int i = 0; i < arr_payload->n_atom_names; i++) {
        const char *name = arr_payload->atom_names[i];
        uint32_t name_len = (uint32_t)strlen(name) + 1; // include null terminator
        status = write_array(CAF_ATOM_NAME_STR, name_len, name, sizeof(char), p_payload,
                             p_payload_bytes);
        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to write atom name string for %s\n", name);
            return;
        }
    }
}

CheckpointStatus validate_array_header(const uint8_t *header, uint16_t *out_flag, uint32_t *out_n)
{
    // TODO: remove array header magic, rely on header flags to indicate what arrays are present
    uint16_t arr_magic;
    memcpy(&arr_magic, header, sizeof(arr_magic));
    if (arr_magic != CHECKPOINT_ARRAY_MAGIC) {
        fprintf(stderr, "Checkpoint error: invalid array magic: expected 0x%X, got 0x%X\n",
                CHECKPOINT_ARRAY_MAGIC, arr_magic);
        return CHECKPOINT_ERROR;
    }
    memcpy(out_flag, header + sizeof(arr_magic), sizeof(uint16_t));
    if (*out_flag <= 0 || *out_flag >= CAF_SENTINEL) {
        fprintf(stderr, "Checkpoint error: invalid array flag: got %u\n", *out_flag);
        return CHECKPOINT_ERROR;
    }
    memcpy(out_n, header + sizeof(arr_magic) + sizeof(uint16_t), sizeof(uint32_t));
    return CHECKPOINT_OK;
}

CheckpointStatus read_array(const uint8_t *payload, size_t *bytes_read, uint16_t *out_flag,
                            uint32_t *out_n, void **out_arr)
{
    CheckpointStatus status = validate_array_header(payload, out_flag, out_n);
    if (status != CHECKPOINT_OK) {
        return status;
    }
    *bytes_read = CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(*out_flag) + sizeof(*out_n);

    size_t elem_size;
    switch (*out_flag) {
    case CAF_ATOMS:
        elem_size = sizeof(CheckpointAtomPayload);
        break;
    case CAF_SUBSTRATE_COMPOSITION:
        elem_size = sizeof(double);
        break;
    case CAF_NN_ENERGY:
        elem_size = sizeof(double);
        break;
    case CAF_IS_SOLUBLE:
        elem_size = sizeof(bool);
        break;
    case CAF_ATOM_NAMES:
        // caller needs to handle this separately since it's an array of strings
        elem_size = sizeof(char *);
        return CHECKPOINT_OK;
        break;
    case CAF_ATOM_NAME_STR:
        elem_size = sizeof(char);
        break;
    default:
        fprintf(stderr, "Checkpoint error: unknown array flag: %u\n", *out_flag);
        return CHECKPOINT_ERROR;
    }

    if (*out_n > 0) {
        *out_arr = malloc((size_t)(*out_n) * elem_size);
        if (*out_arr == NULL) {
            fprintf(stderr, "Checkpoint error: failed to allocate memory for array: %s\n",
                    strerror(errno));
            return CHECKPOINT_ERROR;
        }
        memcpy(*out_arr, payload + *bytes_read, (size_t)(*out_n) * elem_size);
    } else {
        *out_arr = NULL;
    }

    *bytes_read += (size_t)(*out_n) * elem_size;
    // ENHANCE: add some check for validity of the array contents based on the flag?
    return CHECKPOINT_OK;
}

void read_env_arrays(uint8_t *payload, size_t *total_bytes_read,
                     CheckpointEnvArrPayload *arr_payload)
{
    CheckpointStatus status;
    uint16_t flag = 0;
    uint32_t n = 0;
    void *arr = NULL;
    size_t bytes_read = 0;
    uint8_t *payload_ptr = payload;

    // the order of these reads must match the order of writes in fill_env_arrays

    // substrate composition
    status = read_array(payload_ptr, &bytes_read, &flag, &n, &arr);
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
    status = read_array(payload_ptr, &bytes_read, &flag, &n, &arr);
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
    status = read_array(payload_ptr, &bytes_read, &flag, &n, &arr);
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
    status = read_array(payload_ptr, &bytes_read, &flag, &n, &arr);
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
        status = read_array(payload_ptr, &bytes_read, &name_flag, &name_n, (void **)&name_arr);
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

void apply_env_arrays(CheckpointEnvArrPayload *arr_payload, struct SimulationEnv *se)
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
void fill_atom_payload(CheckpointAtomPayload *payload, const Atom *atom)
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
void apply_atom_payload(const CheckpointAtomPayload *payload, Atom *atom)
{
    atom->type = payload->type;
    atom->energy = payload->energy;
    atom->lattice[0] = payload->lattice_u;
    atom->lattice[1] = payload->lattice_v;
    atom->lattice[2] = payload->lattice_w;
    atom->bsradius = payload->bsradius;
}

void write_atom_array(const Atom **atom_arr, uint32_t n_atoms, uint8_t **p_payload,
                      uint32_t *p_payload_bytes)
{
    // write the atom array header first
    write_array((uint16_t)CAF_ATOMS, n_atoms, NULL, 0, p_payload, p_payload_bytes);

    // write each atom's data as a compact payload
    for (uint32_t i = 0; i < n_atoms; ++i) {
        CheckpointAtomPayload atom_payload;
        fill_atom_payload(&atom_payload, atom_arr[i]);
        append_to_payload(&atom_payload, sizeof(CheckpointAtomPayload), p_payload, p_payload_bytes);
    }
}

CheckpointStatus read_atom_array(const uint8_t *payload, size_t *total_bytes_read, uint32_t *out_n,
                                 CheckpointAtomPayload **out_atom_payload_arr)
{
    CheckpointStatus status;
    uint16_t flag = 0;
    uint32_t n = 0;

    status = read_array(payload, total_bytes_read, &flag, &n, (void **)out_atom_payload_arr);
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

CheckpointStatus apply_atom_array(const CheckpointAtomPayload *atom_payload_arr, uint32_t n_atoms,
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

        apply_atom_payload(&atom_payload_arr[i], atom);
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

typedef struct {
    int filler;
} CheckpointOutputFormatArrPayload;

typedef struct {
    int filler;
} CheckpointOutputFormatPayload;

void fill_logging_payload(CheckpointLoggingPayload *payload, const struct LoggingState *ls)
{
    if (!payload || !ls) {
        return;
    }

    payload->framenum = ls->framenum;
    payload->verbose = ls->verbose;
    payload->verbose_interval = ls->verbose_interval;
    payload->increment_precision = ls->increment_precision;
    payload->stime_precision = ls->stime_precision;
    payload->overpot_precision = ls->overpot_precision;
}

CheckpointStatus fill_output_format_array_payload(CheckpointOutputFormatArrPayload *arr_payload,
                                                  const struct LoggingState *ls)
{
    return CHECKPOINT_ERROR;
}

CheckpointStatus write_output_format_array(const CheckpointOutputFormatArrPayload *arr_payload,
                                           uint8_t **p_payload, uint32_t *p_payload_bytes)
{
    return CHECKPOINT_ERROR;
}

CheckpointStatus read_output_format_array(uint8_t *payload, size_t *total_bytes_read,
                                          CheckpointOutputFormatArrPayload *arr_payload)
{
    return CHECKPOINT_ERROR;
}

void apply_logging_payload_to_loggingstate(const CheckpointLoggingPayload *payload,
                                           struct LoggingState *ls)
{
    if (!payload || !ls) {
        return;
    }

    ls->framenum = payload->framenum;
    ls->verbose = payload->verbose;
    ls->verbose_interval = payload->verbose_interval;
    ls->increment_precision = payload->increment_precision;
    ls->stime_precision = payload->stime_precision;
    ls->overpot_precision = payload->overpot_precision;
}

void apply_output_format_array_to_loggingstate(const CheckpointOutputFormatArrPayload *arr_payload,
                                               struct LoggingState *ls)
{
}

/**
 * @brief Create rolling 32-bit checksum for checkpoint file to ensure no corruption
 *
 * @param data pointer to data to checksum
 * @param size size of data in bytes
 * @return uint32_t checksum
 */
uint32_t checkpoint_checksum32(const void *data, size_t size)
{
    const uint8_t *bytes = (const unsigned char *)data;
    uint32_t checksum = 0u;

    // rolling 32-bit checksum - depends on value and order of bytes
    for (size_t i = 0; i < size; ++i) {
        // rotate the checksum by 5 bytes, with bytes pushed out the left being cycled back on
        // the right 5 is ~arbitrary
        checksum = (checksum << 5) | (checksum >> 27);
        // XOR with the next byte
        checksum ^= (uint32_t)bytes[i];
    }

    return checksum;
}
void checkpoint_header_set_magic(CheckpointHeader *header)
{
    memcpy(header->magic, CHECKPOINT_MAGIC_STRING, CHECKPOINT_MAGIC_SIZE);
}

void checkpoint_header_init(CheckpointHeader *header)
{
    memset(header, 0, sizeof(*header));
    checkpoint_header_set_magic(header);
    header->format_version = CHECKPOINT_FORMAT_VERSION;
    header->header_bytes = (uint32_t)sizeof(*header);
}

void checkpoint_header_finalize(CheckpointHeader *header, uint32_t payload_bytes, uint32_t checksum)
{
    header->payload_bytes = payload_bytes;
    header->checksum = checksum;
}

bool checkpoint_header_has_valid_magic(const CheckpointHeader *header)
{
    return memcmp(header->magic, CHECKPOINT_MAGIC_STRING, CHECKPOINT_MAGIC_SIZE) == 0;
}

bool checkpoint_header_has_valid_structure(const CheckpointHeader *header)
{
    return header->format_version == CHECKPOINT_FORMAT_VERSION &&
           header->header_bytes == sizeof(*header) && header->reserved0 == 0u;
}

bool checkpoint_header_is_valid(const CheckpointHeader *header)
{
    if (!checkpoint_header_has_valid_magic(header)) {
        fprintf(stderr, "Checkpoint error: invalid magic in checkpoint header\n");
        return false;
    }
    if (!checkpoint_header_has_valid_structure(header)) {
        fprintf(stderr, "Checkpoint error: invalid shape in checkpoint header\n");
        return false;
    }
    return true;
}

bool checkpoint_header_has_valid_checksum(const CheckpointHeader *header, const void *payload,
                                          size_t payload_size)
{
    return header->checksum == checkpoint_checksum32(payload, payload_size);
}

CheckpointStatus checkpoint_header_write(FILE *file, const CheckpointHeader *header)
{
    size_t written = fwrite(header, sizeof(*header), 1, file);

    if (written != 1) {
        fprintf(stderr, "Checkpoint error: failed to write checkpoint header: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    return CHECKPOINT_OK;
}

CheckpointStatus checkpoint_header_read(FILE *file, CheckpointHeader *header)
{
    size_t read_count = fread(header, sizeof(*header), 1, file);

    if (read_count != 1) {
        fprintf(stderr, "Checkpoint error: failed to read checkpoint header: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    return CHECKPOINT_OK;
}

/**
 * @brief Save a checkpoint file with a header and payload containing the simulation state and
 * environment.
 *
 * @param path file path to write checkpoint to
 * @param ss simulation state to serialize into checkpoint payload (if NULL, only header is
 * written)
 * @param se simulation environment to serialize into checkpoint payload (if NULL, only header
 * is written)
 * @return CheckpointStatus
 */
CheckpointStatus checkpoint_save(const char *path, const struct SimulationState *ss,
                                 const struct SimulationEnv *se, const struct LoggingState *ls)
{
    // create temporary file with the same filename plus a .tmp suffix so interrupted writes do
    // not replace the previous checkpoint.
    char temp_path[512];
    snprintf(temp_path, sizeof(temp_path), "%s.tmp", path);
    FILE *tmp_file = fopen(temp_path, "wb");
    if (!tmp_file) {
        fprintf(stderr, "Checkpoint error: failed to open temp file %s: %s\n", temp_path,
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    // initialize variables
    CheckpointHeader header = {0};
    checkpoint_header_init(&header);

    uint32_t payload_bytes = 0u;
    uint8_t *payload = NULL;
    CheckpointStatus status;

    /* ---- Add primitives to payload ---- */
    CheckpointStatePayload state_payload = {0};
    if (ss) {
        header.has_ss = 1;
        fill_state_payload(&state_payload, ss);
        status = append_to_payload(&state_payload, sizeof(state_payload), &payload, &payload_bytes);
        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to build state payload for %s\n", temp_path);
            fclose(tmp_file);
            remove(temp_path);
            return status;
        }
    }

    CheckpointEnvPayload env_payload = {0};
    if (se) {
        header.has_se = 1;
        fill_env_payload(&env_payload, se);
        status = append_to_payload(&env_payload, sizeof(env_payload), &payload, &payload_bytes);
        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to build state payload for %s\n", temp_path);
            fclose(tmp_file);
            remove(temp_path);
            return status;
        }
    }

    // TODO: logging state
    CheckpointLoggingPayload logging_payload = {0};
    if (ls) {
        header.has_ls = 1;
        fill_logging_payload(&logging_payload, ls);
        status =
            append_to_payload(&logging_payload, sizeof(logging_payload), &payload, &payload_bytes);
        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to build logging payload for %s\n",
                    temp_path);
            fclose(tmp_file);
            remove(temp_path);
            return status;
        }
    }

    /* ---- Add arrays to payload ---- */
    CheckpointEnvArrPayload env_arr_payload = {0};
    if (se) {
        fill_env_array_payload(&env_arr_payload, se);
        write_env_arrays(&env_arr_payload, &payload, &payload_bytes);
        free(env_arr_payload.n_atom_names_str);
    }

    CheckpointOutputFormatArrPayload output_format_arr_payload = {0};
    if (ls) {
        fill_output_format_array_payload(&output_format_arr_payload, ls);
        write_output_format_array(&output_format_arr_payload, &payload, &payload_bytes);
    }

    // if state is present and has atoms, include atom count and array in payload.
    if (ss && ss->atom_cnt > 0) {
        header.has_atoms = 1;
        write_atom_array((const Atom **)ss->atom_arr, (uint32_t)ss->atom_cnt, &payload,
                         &payload_bytes);
    }

    // TODO: add out_formats array to payload

    /* ---- Compute checksum ---- */
    // compute checksum over the combined payload
    uint32_t checksum = checkpoint_checksum32(payload, payload_bytes);
    checkpoint_header_finalize(&header, payload_bytes, checksum);

    /* ---- Write header and payload to temp file ---- */
    // write the fixed-size header first
    status = checkpoint_header_write(tmp_file, &header);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write header to %s\n", temp_path);
        fclose(tmp_file);
        remove(temp_path);
        return status;
    }

    // write the variable-size payload after the header
    if (fwrite(payload, sizeof(*payload), payload_bytes, tmp_file) != payload_bytes) {
        fprintf(stderr, "Checkpoint error: failed to write payload to %s: %s\n", temp_path,
                strerror(errno));
        fclose(tmp_file);
        remove(temp_path);
        return CHECKPOINT_ERROR;
    }

    /* ---- Move temp file to checkpoint path ---- */
    if (fclose(tmp_file) != 0) {
        fprintf(stderr, "Checkpoint error: failed to finalize temp file %s: %s\n", temp_path,
                strerror(errno));
        remove(temp_path);
        return CHECKPOINT_ERROR;
    }

    remove(path);
    if (rename(temp_path, path) != 0) {
        fprintf(stderr, "Checkpoint error: failed to rename %s to %s: %s\n", temp_path, path,
                strerror(errno));
        remove(temp_path);
        return CHECKPOINT_ERROR;
    }

    return CHECKPOINT_OK;
}

// rebuild rate and transition arrays from restored atoms
// must be called after zones are rebuilt
CheckpointStatus rebuild_rates_and_transitions(struct SimulationState *ss, struct SimulationEnv *se)
{
    if (!ss || !ss->atom_arr || ss->atom_cnt <= 0) {
        return CHECKPOINT_OK; // nothing to rebuild
    }

    // reset rate and transition counts to 0 (they will be rebuilt)
    ss->rate_cnt = 0;
    ss->transition_cnt = 0;

    // rebuild transitions for each atom by recomputing its possible moves
    for (long int i = 0; i < ss->atom_cnt; ++i) {
        refresh_transitions(i, ss, se);
    }

    // rebuild the transition probability structure for efficient event selection
    compute_transition_array(ss, se);

    return CHECKPOINT_OK;
}

CheckpointStatus verify_payload_size(FILE *file, CheckpointHeader *header)
{
    if (fseek(file, 0, SEEK_END) != 0) {
        fprintf(stderr, "Checkpoint error: failed to seek to end: %s\n", strerror(errno));
        fclose(file);
        return CHECKPOINT_ERROR;
    }
    long int file_size = ftell(file);
    if (file_size < 0) {
        fprintf(stderr, "Checkpoint error: failed to get size: %s\n", strerror(errno));
        fclose(file);
        return CHECKPOINT_ERROR;
    }
    // confirm payload size matches header
    if (file_size != (long int)(header->header_bytes + header->payload_bytes)) {
        fprintf(stderr, "Checkpoint error: file size mismatch: expected %u bytes, got %ld bytes\n",
                header->header_bytes + header->payload_bytes, file_size);
        fclose(file);
        return CHECKPOINT_ERROR;
    }

    // move file position back to start of payload
    if (fseek(file, (long int)header->header_bytes, SEEK_SET) != 0) {
        fprintf(stderr, "Checkpoint error: failed to seek to payload: %s\n", strerror(errno));
        fclose(file);
        return CHECKPOINT_ERROR;
    }

    return CHECKPOINT_OK;
}

/**
 * @brief Load a checkpoint file, validate its header and checksum, and apply the saved state
 * and environment to the provided SimulationState and SimulationEnv structs.
 *
 * @param path file path to read checkpoint from
 * @param ss simulation state struct to populate from checkpoint payload (if NULL, state payload
 * is ignored)
 * @param se simulation environment struct to populate from checkpoint payload (if NULL, env
 * payload is ignored)
 * @return CheckpointStatus
 */
CheckpointStatus checkpoint_load(const char *path, struct SimulationState *ss,
                                 struct SimulationEnv *se, const struct LoggingState *ls)
{

    // load the checkpoint file
    // read the header
    // validate the header (magic, version, reserved0)
    // get file size and confirm payload size matches header
    // load payload and confirm checksum matches header
    // parse payload in order (state, env, arrays, atoms)
    // apply payload data to SimulationState and SimulationEnv?
    // rebuild derived data and structures (SimEnv arrays, zones, transitions)
    // cleanup and return

    // open the checkpoint and read the header first so we know what payloads follow.
    FILE *file = fopen(path, "rb");
    if (!file) {
        fprintf(stderr, "Checkpoint error: failed to open checkpoint file %s: %s\n", path,
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    CheckpointHeader header;
    CheckpointStatus status = checkpoint_header_read(file, &header);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to read header from %s\n", path);
        fclose(file);
        return status;
    }

    if (!checkpoint_header_is_valid(&header)) {
        fclose(file);
        return CHECKPOINT_ERROR;
    }

    // handle header-only checkpoints (no payload).
    if (header.payload_bytes == 0u) {
        fclose(file);
        return CHECKPOINT_OK;
    }

    // get file size
    status = verify_payload_size(file, &header);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to verify payload size for %s\n", path);
        return status;
    }

    // read in the entire payload
    uint8_t *payload = (uint8_t *)malloc(header.payload_bytes);
    if (!payload) {
        fprintf(stderr, "Checkpoint error: out of memory allocating payload buffer for %s\n", path);
        fclose(file);
        return CHECKPOINT_ERROR;
    }
    size_t read = fread(payload, sizeof(*payload), header.payload_bytes, file);
    if (read != header.payload_bytes) {
        fprintf(stderr, "Checkpoint error: failed to read payload from %s: %s\n", path,
                strerror(errno));
        free(payload);
        fclose(file);
        return CHECKPOINT_ERROR;
    }

    bool checksum_valid = header.checksum == checkpoint_checksum32(payload, header.payload_bytes);
    if (!checksum_valid) {
        fprintf(stderr, "Checkpoint error: checksum mismatch in %s: expected %u, got %u\n", path,
                header.checksum, checkpoint_checksum32(payload, header.payload_bytes));
        free(payload);
        fclose(file);
        return CHECKPOINT_ERROR;
    }

    // payload is valid, we are done with the file
    fclose(file);

    // parse payload in order (state, env, arrays, atoms)
    uint8_t *payload_ptr = payload;

    // state payload
    CheckpointStatePayload state_payload = {0};
    if ((header.has_ss || header.has_atoms) && !ss) {
        fprintf(stderr,
                "Checkpoint error: no SimulationState object provided to load state payload from "
                "checkpoint %s\n",
                path);
        free(payload);
        return CHECKPOINT_ERROR;
    }
    if (header.has_ss) {
        if (payload_ptr + sizeof(CheckpointStatePayload) > payload + header.payload_bytes) {
            fprintf(stderr, "Checkpoint error: payload too small for state payload in %s\n", path);
            free(payload);
            return CHECKPOINT_ERROR;
        }
        memcpy(&state_payload, payload_ptr, sizeof(CheckpointStatePayload));
        payload_ptr += sizeof(CheckpointStatePayload);
    }

    // env payload
    CheckpointEnvPayload env_payload = {0};
    if (header.has_se && !se) {
        fprintf(stderr,
                "Checkpoint error: no SimulationEnv object provided to load env payload from "
                "checkpoint %s\n",
                path);
        free(payload);
        return CHECKPOINT_ERROR;
    }
    if (header.has_se) {
        if (payload_ptr + sizeof(CheckpointEnvPayload) > payload + header.payload_bytes) {
            fprintf(stderr, "Checkpoint error: payload too small for env payload in %s\n", path);
            free(payload);
            return CHECKPOINT_ERROR;
        }
        memcpy(&env_payload, payload_ptr, sizeof(CheckpointEnvPayload));
        payload_ptr += sizeof(CheckpointEnvPayload);
    }

    // logging payload
    CheckpointLoggingPayload logging_payload = {0};
    if (header.has_ls && !ls) {
        fprintf(stderr,
                "Checkpoint error: no SimulationEnv object provided to load env payload from "
                "checkpoint %s\n",
                path);
        free(payload);
        return CHECKPOINT_ERROR;
    }
    if (header.has_ls) {
        if (payload_ptr + sizeof(CheckpointLoggingPayload) > payload + header.payload_bytes) {
            fprintf(stderr, "Checkpoint error: payload too small for env payload in %s\n", path);
            free(payload);
            return CHECKPOINT_ERROR;
        }
        memcpy(&logging_payload, payload_ptr, sizeof(CheckpointLoggingPayload));
        payload_ptr += sizeof(CheckpointLoggingPayload);
    }

    // SimEnv arrays
    CheckpointEnvArrPayload env_arr_payload = {0};
    if (header.has_se) {
        size_t bytes_read;
        read_env_arrays(payload_ptr, &bytes_read, &env_arr_payload);
        payload_ptr += bytes_read;
    }

    // output format arrays
    CheckpointOutputFormatArrPayload output_format_arr_payload = {0};
    if (header.has_ls) {
        size_t bytes_read;
        read_output_format_array(payload_ptr, &bytes_read, &output_format_arr_payload);
        payload_ptr += bytes_read;
    }

    // Atom array
    CheckpointAtomPayload *atom_payload_arr;
    uint32_t atom_cnt;
    if (header.has_atoms) {
        size_t bytes_read;
        read_atom_array(payload_ptr, &bytes_read, &atom_cnt, &atom_payload_arr);
        payload_ptr += bytes_read;
    }

    // apply payloads
    // apply ss and se first
    // calculate derived members
    // apply se arrays
    // initialize simulation arrays
    // populate atom and zone arrays
    // initialize transitions and rates

    // cases where objects are required but not provided should have been caught above
    struct SimulationConfig config = {0};
    if (header.has_ss) {
        apply_state_payload_to_simstate(&state_payload, ss);
    }
    if (header.has_se) {
        apply_env_payload_to_config(&env_payload, &config);
        initialize_env_from_config(&config, se);

        apply_env_arrays(&env_arr_payload, se);
    }
    // TODO: apply logging payload and output format arrays
    if (header.has_ss && header.has_se) {
        allocate_simulation_arrays(ss, se);
        if (header.has_atoms) {
            errno = 0;
            status = apply_atom_array(atom_payload_arr, atom_cnt, env_payload.num_elements, ss, se);
            if (status != CHECKPOINT_OK) {
                fprintf(stderr,
                        "Checkpoint error: failed to apply atom payload to atom array: %s\n",
                        strerror(errno));
                // TODO: free malloc'd arrays on failure
                return status;
            }
        }
        rebuild_rates_and_transitions(ss, se);
    }

    // env arrays in env_arr_payload have ownership transferred to SimulationEnv
    if (header.has_atoms && atom_cnt > 0) {
        free(atom_payload_arr);
    }

    free(payload);
    return CHECKPOINT_OK;
}
