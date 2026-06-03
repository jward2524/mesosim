#include "Checkpoint.h"
#include "CheckpointLogging.h"
#include "CheckpointSimulation.h"
#include "CheckpointUtils.h"
#include "Initialization.h"
#include "Simulation.h"
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
void set_checkpoint_header_magic(CheckpointHeader *header)
{
    memcpy(header->magic, CHECKPOINT_MAGIC_STRING, CHECKPOINT_MAGIC_SIZE);
}

void initialize_checkpoint_header(CheckpointHeader *header)
{
    memset(header, 0, sizeof(*header));
    set_checkpoint_header_magic(header);
    header->format_version = CHECKPOINT_FORMAT_VERSION;
    header->header_bytes = (uint32_t)sizeof(*header);
}

void finalize_checkpoint_header(CheckpointHeader *header, uint32_t payload_bytes, uint32_t checksum)
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

CheckpointStatus write_checkpoint_header(FILE *file, const CheckpointHeader *header)
{
    size_t written = fwrite(header, sizeof(*header), 1, file);

    if (written != 1) {
        fprintf(stderr, "Checkpoint error: failed to write checkpoint header: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    return CHECKPOINT_OK;
}

CheckpointStatus write_checkpoint_payload(FILE *file, const uint8_t *payload,
                                          const uint32_t payload_bytes)
{
    size_t written = fwrite(payload, sizeof(*payload), payload_bytes, file);

    if (written != payload_bytes) {
        fprintf(stderr, "Checkpoint error: failed to write checkpoint header: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    return CHECKPOINT_OK;
}

CheckpointStatus read_checkpoint_header(FILE *file, CheckpointHeader *header)
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
CheckpointStatus write_checkpoint_buffer(const char *path, const struct SimulationState *ss,
                                         const struct SimulationEnv *se,
                                         const struct LoggingState *ls)
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
    initialize_checkpoint_header(&header);

    uint32_t payload_bytes = 0u;
    uint8_t *payload = NULL;
    CheckpointStatus status;

    /* ---- Add primitives to payload ---- */
    SimStatePayload state_payload = {0};
    if (ss) {
        header.has_ss = 1;
        pack_simstate(&state_payload, ss);
        status = append_to_payload(&state_payload, sizeof(state_payload), &payload, &payload_bytes);
        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to build state payload for %s\n", temp_path);
            fclose(tmp_file);
            remove(temp_path);
            return status;
        }
    }

    SimEnvPayload env_payload = {0};
    if (se) {
        header.has_se = 1;
        pack_simenv(&env_payload, se);
        status = append_to_payload(&env_payload, sizeof(env_payload), &payload, &payload_bytes);
        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to build state payload for %s\n", temp_path);
            fclose(tmp_file);
            remove(temp_path);
            return status;
        }
    }

    // TODO: logging state
    LoggingPayload logging_payload = {0};
    if (ls) {
        header.has_ls = 1;
        pack_logging_payload(&logging_payload, ls);
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
    SimEnvArrPayload env_arr_payload = {0};
    if (se) {
        pack_env_arrays(&env_arr_payload, se);
        serialize_env_arrays(&env_arr_payload, &payload, &payload_bytes);
        free(env_arr_payload.n_atom_names_str);
    }

    // write output formats array even if it is empty
    OutFormatArrPayload output_format_arr_payload = {0};
    if (ls) {
        pack_output_format_array(&output_format_arr_payload, ls);
        serialize_output_format_array(&output_format_arr_payload, &payload, &payload_bytes);
        free(output_format_arr_payload.formats);
    }

    // if state is present and has atoms, include atom count and array in payload.
    if (ss && ss->atom_cnt > 0) {
        header.has_atoms = 1;
        serialize_atom_array((const Atom **)ss->atom_arr, (uint32_t)ss->atom_cnt, &payload,
                             &payload_bytes);
    }

    /* ---- Compute checksum ---- */
    // compute checksum over the combined payload
    uint32_t checksum = checkpoint_checksum32(payload, payload_bytes);
    finalize_checkpoint_header(&header, payload_bytes, checksum);

    /* ---- Write header and payload to temp file ---- */
    // write the fixed-size header first
    status = write_checkpoint_header(tmp_file, &header);
    if (status != CHECKPOINT_OK) {
        fclose(tmp_file);
        remove(temp_path);
        return status;
    }

    // write the variable-size payload after the header
    status = write_checkpoint_payload(tmp_file, payload, payload_bytes);
    if (status != CHECKPOINT_OK) {
        fclose(tmp_file);
        remove(temp_path);
        return status;
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
CheckpointStatus read_checkpoint_file(const char *path, struct SimulationState *ss,
                                      struct SimulationEnv *se, struct LoggingState *ls)
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
    CheckpointStatus status = read_checkpoint_header(file, &header);
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

    // TODO: make an unserialize_checkpoint_payload function
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
    SimStatePayload state_payload = {0};
    if ((header.has_ss || header.has_atoms) && !ss) {
        fprintf(stderr,
                "Checkpoint error: no SimulationState object provided to load state payload from "
                "checkpoint %s\n",
                path);
        free(payload);
        return CHECKPOINT_ERROR;
    }
    if (header.has_ss) {
        if (payload_ptr + sizeof(SimStatePayload) > payload + header.payload_bytes) {
            fprintf(stderr, "Checkpoint error: payload too small for state payload in %s\n", path);
            free(payload);
            return CHECKPOINT_ERROR;
        }
        memcpy(&state_payload, payload_ptr, sizeof(SimStatePayload));
        payload_ptr += sizeof(SimStatePayload);
    }

    // env payload
    SimEnvPayload env_payload = {0};
    if (header.has_se && !se) {
        fprintf(stderr,
                "Checkpoint error: no SimulationEnv object provided to load env payload from "
                "checkpoint %s\n",
                path);
        free(payload);
        return CHECKPOINT_ERROR;
    }
    if (header.has_se) {
        if (payload_ptr + sizeof(SimEnvPayload) > payload + header.payload_bytes) {
            fprintf(stderr, "Checkpoint error: payload too small for env payload in %s\n", path);
            free(payload);
            return CHECKPOINT_ERROR;
        }
        memcpy(&env_payload, payload_ptr, sizeof(SimEnvPayload));
        payload_ptr += sizeof(SimEnvPayload);
    }

    // logging payload
    LoggingPayload logging_payload = {0};
    if (header.has_ls && !ls) {
        fprintf(stderr,
                "Checkpoint error: no SimulationEnv object provided to load env payload from "
                "checkpoint %s\n",
                path);
        free(payload);
        return CHECKPOINT_ERROR;
    }
    if (header.has_ls) {
        if (payload_ptr + sizeof(LoggingPayload) > payload + header.payload_bytes) {
            fprintf(stderr, "Checkpoint error: payload too small for env payload in %s\n", path);
            free(payload);
            return CHECKPOINT_ERROR;
        }
        memcpy(&logging_payload, payload_ptr, sizeof(LoggingPayload));
        payload_ptr += sizeof(LoggingPayload);
    }

    // SimEnv arrays
    SimEnvArrPayload env_arr_payload = {0};
    if (header.has_se) {
        size_t bytes_read;
        unserialize_env_arrays(payload_ptr, &bytes_read, &env_arr_payload);
        payload_ptr += bytes_read;
    }

    // output format arrays
    OutFormatArrPayload output_format_arr_payload = {0};
    if (header.has_ls) {
        size_t bytes_read;
        unserialize_output_format_array(payload_ptr, &bytes_read, &output_format_arr_payload);
        payload_ptr += bytes_read;
    }

    // Atom array
    AtomPayload *atom_payload_arr;
    uint32_t atom_cnt;
    if (header.has_atoms) {
        size_t bytes_read;
        unserialize_atom_array(payload_ptr, &bytes_read, &atom_cnt, &atom_payload_arr);
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
        unpack_state_payload_to_simstate(&state_payload, ss);
    }
    if (header.has_se) {
        unpack_simenv_payload_to_config(&env_payload, &config);
        initialize_env_from_config(&config, se);

        unpack_env_arrays(&env_arr_payload, se);
    }

    if (header.has_ls) {
        unpack_logging_payload(&logging_payload, ls);
        unpack_output_format_array(&output_format_arr_payload, ls);
    }

    if (header.has_ss && header.has_se) {
        allocate_simulation_arrays(ss, se);
        if (header.has_atoms) {
            errno = 0;
            status =
                unpack_atom_array(atom_payload_arr, atom_cnt, env_payload.num_elements, ss, se);
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
