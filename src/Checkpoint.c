#include "Checkpoint.h"
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// packed payload used for the first real checkpoint slice: the mutable state scalars that need
// to survive a pause/resume cycle on the same machine.
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
} CheckpointStatePayload;
#pragma pack(pop)

// copy the SimulationState fields we care about into a compact on-disk payload.
static void fill_state_payload(CheckpointStatePayload *payload, const struct SimulationState *ss)
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
static void apply_state_payload(const CheckpointStatePayload *payload, struct SimulationState *ss)
{
    ss->iter = payload->iter;
    ss->mmc_steps = payload->mmc_steps;
    ss->final_iteration = payload->final_iteration;
    ss->run_stime = payload->run_stime;
    ss->simulation_should_kill_itself = payload->simulation_should_kill_itself;
    ss->elapsed_stime = payload->elapsed_stime;
    ss->sim_end_type = payload->sim_end_type;
    ss->frequency_sum = payload->frequency_sum;
    ss->total_internal_energy = payload->total_internal_energy;
    ss->temperature = payload->temperature;
    ss->overpotential = payload->overpotential;
    ss->total_atoms_dissolved = payload->total_atoms_dissolved;
}

// packed payload for the immutable simulation environment: geometry, lattice, and configuration.
// these fields define how to interpret the atom array and rebuild derived structures on restore.
#pragma pack(push, 1)
typedef struct {
    unsigned flavor;
    unsigned rand_seed;
    int geometry;
    bool evaporation_flag;
    int lattice_type;
    int sheet_thickness;
    int cluster_radius;
    int system_size_x;
    int system_size_y;
    int system_size_z;
    int num_elements;
    int num_nn_levels;
    int num_bond_types;
    int num_neighbor_types;
    int num_transition_vectors;
    double primitive_basis_u1;
    double primitive_basis_u2;
    double primitive_basis_u3;
    double primitive_basis_v1;
    double primitive_basis_v2;
    double primitive_basis_v3;
    double primitive_basis_w1;
    double primitive_basis_w2;
    double primitive_basis_w3;
    double initial_overpotential;
    double overpotential_ramp_rate;
    double max_overpotential;
} CheckpointEnvPayload;
#pragma pack(pop)

// copy the immutable SimulationEnv fields into a compact on-disk payload for restore validation.
static void fill_env_payload(CheckpointEnvPayload *payload, const struct SimulationEnv *se)
{
    payload->flavor = se->flavor;
    payload->rand_seed = se->rand_seed;
    payload->geometry = se->geometry;
    payload->evaporation_flag = se->evaporation_flag;
    payload->lattice_type = se->lattice_type;
    payload->sheet_thickness = se->sheet_thickness;
    payload->cluster_radius = se->cluster_radius;
    payload->system_size_x = se->system_size_x;
    payload->system_size_y = se->system_size_y;
    payload->system_size_z = se->system_size_z;
    payload->num_elements = se->num_elements;
    payload->num_nn_levels = se->num_nn_levels;
    payload->num_bond_types = se->num_bond_types;
    payload->num_neighbor_types = se->num_neighbor_types;
    payload->num_transition_vectors = se->num_transition_vectors;
    // copy 3x3 primitive basis matrix
    payload->primitive_basis_u1 = se->primitive_basis[0][0];
    payload->primitive_basis_u2 = se->primitive_basis[0][1];
    payload->primitive_basis_u3 = se->primitive_basis[0][2];
    payload->primitive_basis_v1 = se->primitive_basis[1][0];
    payload->primitive_basis_v2 = se->primitive_basis[1][1];
    payload->primitive_basis_v3 = se->primitive_basis[1][2];
    payload->primitive_basis_w1 = se->primitive_basis[2][0];
    payload->primitive_basis_w2 = se->primitive_basis[2][1];
    payload->primitive_basis_w3 = se->primitive_basis[2][2];
    payload->initial_overpotential = se->initial_overpotential;
    payload->overpotential_ramp_rate = se->overpotential_ramp_rate;
    payload->max_overpotential = se->max_overpotential;
}

// restore the immutable SimulationEnv config from a saved payload.
static void apply_env_payload(const CheckpointEnvPayload *payload, struct SimulationEnv *se)
{
    se->flavor = payload->flavor;
    se->rand_seed = payload->rand_seed;
    se->geometry = payload->geometry;
    se->evaporation_flag = payload->evaporation_flag;
    se->lattice_type = payload->lattice_type;
    se->sheet_thickness = payload->sheet_thickness;
    se->cluster_radius = payload->cluster_radius;
    se->system_size_x = payload->system_size_x;
    se->system_size_y = payload->system_size_y;
    se->system_size_z = payload->system_size_z;
    se->num_elements = payload->num_elements;
    se->num_nn_levels = payload->num_nn_levels;
    se->num_bond_types = payload->num_bond_types;
    se->num_neighbor_types = payload->num_neighbor_types;
    se->num_transition_vectors = payload->num_transition_vectors;
    // restore 3x3 primitive basis matrix
    se->primitive_basis[0][0] = payload->primitive_basis_u1;
    se->primitive_basis[0][1] = payload->primitive_basis_u2;
    se->primitive_basis[0][2] = payload->primitive_basis_u3;
    se->primitive_basis[1][0] = payload->primitive_basis_v1;
    se->primitive_basis[1][1] = payload->primitive_basis_v2;
    se->primitive_basis[1][2] = payload->primitive_basis_v3;
    se->primitive_basis[2][0] = payload->primitive_basis_w1;
    se->primitive_basis[2][1] = payload->primitive_basis_w2;
    se->primitive_basis[2][2] = payload->primitive_basis_w3;
    se->initial_overpotential = payload->initial_overpotential;
    se->overpotential_ramp_rate = payload->overpotential_ramp_rate;
    se->max_overpotential = payload->max_overpotential;
}

// packed payload for a single atom: just enough to recreate it with its position and type.
// derived fields like transition_indices, neighbor_atom_idxs, and linked-list pointers will be
// rebuilt from the atom snapshot during restore, so we do not serialize them here.
#pragma pack(push, 1)
typedef struct {
    unsigned char type;
    double energy;
    int lattice_u;
    int lattice_v;
    int lattice_w;
    double cartesian_x;
    double cartesian_y;
    double cartesian_z;
    double bsradius;
} CheckpointAtomPayload;
#pragma pack(pop)

// serialize a single atom into its checkpoint representation.
static void fill_atom_payload(CheckpointAtomPayload *payload, const Atom *atom)
{
    payload->type = atom->type;
    payload->energy = atom->energy;
    payload->lattice_u = atom->lattice[0];
    payload->lattice_v = atom->lattice[1];
    payload->lattice_w = atom->lattice[2];
    payload->cartesian_x = atom->cartesian[0];
    payload->cartesian_y = atom->cartesian[1];
    payload->cartesian_z = atom->cartesian[2];
    payload->bsradius = atom->bsradius;
}

// restore a single atom from its checkpoint representation.
// note: transition_indices, neighbor_atom_idxs, and linked-list fields are zeroed;
// they will be rebuilt from the simulation state during restore.
static void apply_atom_payload(const CheckpointAtomPayload *payload, Atom *atom)
{
    atom->type = payload->type;
    atom->energy = payload->energy;
    atom->lattice[0] = payload->lattice_u;
    atom->lattice[1] = payload->lattice_v;
    atom->lattice[2] = payload->lattice_w;
    atom->cartesian[0] = payload->cartesian_x;
    atom->cartesian[1] = payload->cartesian_y;
    atom->cartesian[2] = payload->cartesian_z;
    atom->bsradius = payload->bsradius;
}

/**
 * @brief Create checksum for checkpoint file to ensure no corruption
 *
 * @param data pointer to data to checksum
 * @param size size of data in bytes
 * @return uint32_t checksum
 */
uint32_t checkpoint_checksum32(const void *data, size_t size)
{
    const unsigned char *bytes = (const unsigned char *)data;
    uint32_t checksum = 0u;

    // checksum that depends on value and order of bytes
    for (size_t i = 0; i < size; ++i) {
        // rotate the checksum by 5 bytes, with bytes pushed out the left being cycled back on the
        // right
        // 5 is ~arbitrary
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

bool checkpoint_header_has_valid_shape(const CheckpointHeader *header)
{
    return header->format_version == CHECKPOINT_FORMAT_VERSION &&
           header->header_bytes == sizeof(*header) && header->reserved0 == 0u;
}

bool checkpoint_header_has_valid_checksum(const CheckpointHeader *header, const void *payload,
                                          size_t payload_size)
{
    return header->checksum == checkpoint_checksum32(payload, payload_size);
}

/**
 * @brief Save a checkpoint file with a header and payload containing the simulation state and
 * environment.
 *
 * @param path file path to write checkpoint to
 * @param ss simulation state to serialize into checkpoint payload (if NULL, only header is written)
 * @param se simulation environment to serialize into checkpoint payload (if NULL, only header is
 * written)
 * @return CheckpointStatus
 */
CheckpointStatus checkpoint_save(const char *path, const struct SimulationState *ss,
                                 const struct SimulationEnv *se)
{
    // create temporary file with the same filename plus a .tmp suffix so interrupted writes do not
    // replace the previous checkpoint.
    char temp_path[512];
    snprintf(temp_path, sizeof(temp_path), "%s.tmp", path);
    FILE *file = fopen(temp_path, "wb");
    if (!file) {
        fprintf(stderr, "Checkpoint error: failed to open temp file %s: %s\n", temp_path,
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    // prepare both payloads if available, then calculate combined checksum for file integrity.
    CheckpointHeader header;
    checkpoint_header_init(&header);

    CheckpointStatePayload state_payload;
    CheckpointEnvPayload env_payload;
    uint32_t payload_bytes = 0u;
    uint32_t checksum = 0u;

    // calculate total payload size and combined checksum over both payloads if present.
    if (ss) {
        fill_state_payload(&state_payload, ss);
        payload_bytes += (uint32_t)sizeof(state_payload);
    }
    if (se) {
        fill_env_payload(&env_payload, se);
        payload_bytes += (uint32_t)sizeof(env_payload);
    }

    // if state is present and has atoms, include atom count and array in payload.
    uint32_t atom_payload_bytes = 0u;
    if (ss && ss->atom_cnt > 0) {
        atom_payload_bytes =
            (uint32_t)(sizeof(uint32_t) + ss->atom_cnt * sizeof(CheckpointAtomPayload));
        payload_bytes += atom_payload_bytes;
    }

    // compute checksum over the combined payload: state first, then env, then atoms.
    if (ss) {
        checksum = checkpoint_checksum32(&state_payload, sizeof(state_payload));
    }
    if (se) {
        uint32_t env_checksum = checkpoint_checksum32(&env_payload, sizeof(env_payload));
        if (ss) {
            checksum ^= env_checksum;
        } else {
            checksum = env_checksum;
        }
    }

    // compute checksum for atoms if present.
    if (atom_payload_bytes > 0) {
        uint32_t atom_cnt_val = (uint32_t)ss->atom_cnt;
        uint32_t atom_count_checksum = checkpoint_checksum32(&atom_cnt_val, sizeof(atom_cnt_val));
        for (long int i = 0; i < ss->atom_cnt; ++i) {
            CheckpointAtomPayload atom_payload;
            fill_atom_payload(&atom_payload, ss->atom_arr[i]);
            uint32_t atom_checksum = checkpoint_checksum32(&atom_payload, sizeof(atom_payload));
            atom_count_checksum ^= atom_checksum;
        }
        checksum ^= atom_count_checksum;
    }

    checkpoint_header_finalize(&header, payload_bytes, checksum);

    // write the fixed-size header first.
    CheckpointStatus status = checkpoint_header_write(file, &header);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to write header to %s\n", temp_path);
        fclose(file);
        remove(temp_path);
        return status;
    }

    // append state payload if available.
    if (ss) {
        size_t written = fwrite(&state_payload, sizeof(state_payload), 1, file);
        if (written != 1) {
            fprintf(stderr, "Checkpoint error: failed to write state payload to %s\n", temp_path);
            fclose(file);
            remove(temp_path);
            return CHECKPOINT_ERROR;
        }
    }

    // append env payload if available.
    if (se) {
        size_t written = fwrite(&env_payload, sizeof(env_payload), 1, file);
        if (written != 1) {
            fprintf(stderr, "Checkpoint error: failed to write env payload to %s\n", temp_path);
            fclose(file);
            remove(temp_path);
            return CHECKPOINT_ERROR;
        }
    }

    // append atom array if present: write count first, then each atom.
    if (atom_payload_bytes > 0) {
        uint32_t atom_cnt_val = (uint32_t)ss->atom_cnt;
        size_t written = fwrite(&atom_cnt_val, sizeof(atom_cnt_val), 1, file);
        if (written != 1) {
            fprintf(stderr, "Checkpoint error: failed to write atom count to %s\n", temp_path);
            fclose(file);
            remove(temp_path);
            return CHECKPOINT_ERROR;
        }
        for (long int i = 0; i < ss->atom_cnt; ++i) {
            CheckpointAtomPayload atom_payload;
            fill_atom_payload(&atom_payload, ss->atom_arr[i]);
            written = fwrite(&atom_payload, sizeof(atom_payload), 1, file);
            if (written != 1) {
                fprintf(stderr, "Checkpoint error: failed to write atom %ld to %s\n", i, temp_path);
                fclose(file);
                remove(temp_path);
                return CHECKPOINT_ERROR;
            }
        }
    }

    // move the completed temp file into place.
    if (fclose(file) != 0) {
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

/**
 * @brief Load a checkpoint file, validate its header and checksum, and apply the saved state and
 * environment to the provided SimulationState and SimulationEnv structs.
 *
 * @param path file path to read checkpoint from
 * @param ss simulation state struct to populate from checkpoint payload (if NULL, state payload is
 * ignored)
 * @param se simulation environment struct to populate from checkpoint payload (if NULL, env payload
 * is ignored)
 * @return CheckpointStatus
 */
CheckpointStatus checkpoint_load(const char *path, struct SimulationState *ss,
                                 struct SimulationEnv *se)
{
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

    if (!checkpoint_header_has_valid_magic(&header) ||
        !checkpoint_header_has_valid_shape(&header)) {
        fprintf(stderr, "Checkpoint error: invalid header (magic/shape) in %s\n", path);
        fclose(file);
        return CHECKPOINT_ERROR;
    }

    // handle header-only checkpoints (no payload).
    if (header.payload_bytes == 0u) {
        fclose(file);
        return CHECKPOINT_OK;
    }

    // determine what payloads are present based on total payload_bytes.
    uint32_t state_payload_size = (uint32_t)sizeof(CheckpointStatePayload);
    uint32_t env_payload_size = (uint32_t)sizeof(CheckpointEnvPayload);
    bool has_state = false;
    bool has_env = false;
    bool has_atoms = false;

    // infer which payloads are present from their combined size.
    uint32_t expected_size_without_atoms = 0u;
    if (ss) {
        expected_size_without_atoms += state_payload_size;
        has_state = true;
    }
    if (se) {
        expected_size_without_atoms += env_payload_size;
        has_env = true;
    }

    // if payload_bytes is larger than state+env, atoms are present (only if state is provided).
    if (ss && header.payload_bytes > expected_size_without_atoms) {
        has_atoms = true;
    } else if (header.payload_bytes != expected_size_without_atoms &&
               expected_size_without_atoms > 0u) {
        // size mismatch: can't determine structure.
        fprintf(stderr, "Checkpoint error: unexpected payload size in %s\n", path);
        fclose(file);
        return CHECKPOINT_ERROR;
    } else if (expected_size_without_atoms == 0u) {
        // no payloads expected but payload_bytes > 0: invalid.
        fclose(file);
        return CHECKPOINT_ERROR;
    }

    // read state payload if present.
    CheckpointStatePayload state_payload;
    if (has_state) {
        if (fread(&state_payload, sizeof(state_payload), 1, file) != 1) {
            fprintf(stderr, "Checkpoint error: failed to read state payload from %s\n", path);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
    }

    // read env payload if present.
    CheckpointEnvPayload env_payload;
    if (has_env) {
        if (fread(&env_payload, sizeof(env_payload), 1, file) != 1) {
            fprintf(stderr, "Checkpoint error: failed to read env payload from %s\n", path);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
    }

    // read atom count and array if present.
    uint32_t atom_cnt_loaded = 0u;
    CheckpointAtomPayload *atom_payloads = NULL;
    if (has_atoms) {
        // read atom count first.
        if (fread(&atom_cnt_loaded, sizeof(atom_cnt_loaded), 1, file) != 1) {
            fprintf(stderr, "Checkpoint error: failed to read atom count from %s\n", path);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
        // allocate temp buffer for atom payloads.
        atom_payloads =
            (CheckpointAtomPayload *)malloc(atom_cnt_loaded * sizeof(CheckpointAtomPayload));
        if (!atom_payloads) {
            fprintf(stderr, "Checkpoint error: out of memory allocating atom payload buffer (%u)\n",
                    atom_cnt_loaded);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
        // read all atoms.
        for (uint32_t i = 0u; i < atom_cnt_loaded; ++i) {
            if (fread(&atom_payloads[i], sizeof(CheckpointAtomPayload), 1, file) != 1) {
                fprintf(stderr, "Checkpoint error: failed to read atom %u from %s\n", i, path);
                free(atom_payloads);
                fclose(file);
                return CHECKPOINT_ERROR;
            }
        }
    }

    // verify combined checksum before applying any data.
    if (has_state && has_env && has_atoms) {
        // verify all three together.
        uint32_t state_check = checkpoint_checksum32(&state_payload, sizeof(state_payload));
        uint32_t env_check = checkpoint_checksum32(&env_payload, sizeof(env_payload));
        uint32_t atom_count_check =
            checkpoint_checksum32(&atom_cnt_loaded, sizeof(atom_cnt_loaded));
        for (uint32_t i = 0u; i < atom_cnt_loaded; ++i) {
            uint32_t atom_check =
                checkpoint_checksum32(&atom_payloads[i], sizeof(atom_payloads[i]));
            atom_count_check ^= atom_check;
        }
        if (header.checksum != (state_check ^ env_check ^ atom_count_check)) {
            fprintf(stderr, "Checkpoint error: checksum mismatch in %s (state+env+atoms)\n", path);
            free(atom_payloads);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
    } else if (has_state && has_atoms) {
        // verify state and atoms together.
        uint32_t state_check = checkpoint_checksum32(&state_payload, sizeof(state_payload));
        uint32_t atom_count_check =
            checkpoint_checksum32(&atom_cnt_loaded, sizeof(atom_cnt_loaded));
        for (uint32_t i = 0u; i < atom_cnt_loaded; ++i) {
            uint32_t atom_check =
                checkpoint_checksum32(&atom_payloads[i], sizeof(atom_payloads[i]));
            atom_count_check ^= atom_check;
        }
        if (header.checksum != (state_check ^ atom_count_check)) {
            fprintf(stderr, "Checkpoint error: checksum mismatch in %s (state+atoms)\n", path);
            free(atom_payloads);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
    } else if (has_state && has_env) {
        // verify state and env together.
        uint32_t state_check = checkpoint_checksum32(&state_payload, sizeof(state_payload));
        uint32_t env_check = checkpoint_checksum32(&env_payload, sizeof(env_payload));
        if (header.checksum != (state_check ^ env_check)) {
            fprintf(stderr, "Checkpoint error: checksum mismatch in %s (state+env)\n", path);
            free(atom_payloads);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
    } else if (has_state) {
        // verify state-only checksum.
        if (!checkpoint_header_has_valid_checksum(&header, &state_payload, sizeof(state_payload))) {
            fprintf(stderr, "Checkpoint error: checksum mismatch in %s (state only)\n", path);
            free(atom_payloads);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
    } else if (has_env) {
        // verify env-only checksum.
        if (!checkpoint_header_has_valid_checksum(&header, &env_payload, sizeof(env_payload))) {
            fprintf(stderr, "Checkpoint error: checksum mismatch in %s (env only)\n", path);
            free(atom_payloads);
            fclose(file);
            return CHECKPOINT_ERROR;
        }
    }

    // apply payloads to the provided structures.
    if (has_state && ss) {
        apply_state_payload(&state_payload, ss);
        ss->atom_cnt = atom_cnt_loaded;
        // allocate atom_arr if atoms were loaded.
        if (has_atoms && atom_cnt_loaded > 0u) {
            ss->atom_arr = (Atom **)malloc(atom_cnt_loaded * sizeof(Atom *));
            if (!ss->atom_arr) {
                fprintf(stderr, "Checkpoint error: out of memory allocating atom_arr (%u)\n",
                        atom_cnt_loaded);
                free(atom_payloads);
                fclose(file);
                return CHECKPOINT_ERROR;
            }
            for (long int i = 0; i < (long int)atom_cnt_loaded; ++i) {
                ss->atom_arr[i] = (Atom *)malloc(sizeof(Atom));
                if (!ss->atom_arr[i]) {
                    fprintf(stderr, "Checkpoint error: out of memory allocating atom %ld\n", i);
                    free(atom_payloads);
                    fclose(file);
                    return CHECKPOINT_ERROR;
                }
                // zero-initialize the atom, then apply payload (sets core fields).
                memset(ss->atom_arr[i], 0, sizeof(Atom));
                apply_atom_payload(&atom_payloads[i], ss->atom_arr[i]);
            }
        }
    }
    if (has_env && se) {
        apply_env_payload(&env_payload, se);
    }

    free(atom_payloads);
    fclose(file);
    return CHECKPOINT_OK;
}

CheckpointStatus checkpoint_header_write(FILE *file, const CheckpointHeader *header)
{
    size_t written = fwrite(header, sizeof(*header), 1, file);

    if (written != 1) {
        return CHECKPOINT_ERROR;
    }

    return CHECKPOINT_OK;
}

CheckpointStatus checkpoint_header_read(FILE *file, CheckpointHeader *header)
{
    size_t read_count = fread(header, sizeof(*header), 1, file);

    if (read_count != 1) {
        return CHECKPOINT_ERROR;
    }

    return CHECKPOINT_OK;
}
