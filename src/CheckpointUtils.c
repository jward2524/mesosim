#include "CheckpointUtils.h"
#include "CheckpointLogging.h"
#include "CheckpointSimulation.h"
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

CheckpointStatus serialize_array_magic(uint8_t **p_payload, uint32_t *p_payload_bytes)
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
CheckpointStatus serialize_array_header_into_payload(uint16_t flag, uint32_t n, uint8_t **p_payload)
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
CheckpointStatus serialize_array(uint16_t flag, uint32_t n, const void *arr, size_t elem_size,
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

    status = serialize_array_header_into_payload(flag, n, &intermediate_payload);
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

static bool array_flag_elem_size(uint16_t flag, size_t *elem_size)
{
    if (elem_size == NULL) {
        return false;
    }

    switch (flag) {
    case CAF_ATOMS:
        *elem_size = sizeof(AtomPayload);
        return true;
    case CAF_SUBSTRATE_COMPOSITION:
    case CAF_NN_ENERGY:
        *elem_size = sizeof(double);
        return true;
    case CAF_IS_SOLUBLE:
        *elem_size = sizeof(bool);
        return true;
    case CAF_OUTPUT_FORMATS:
        *elem_size = sizeof(OutFormatPayload);
        return true;
    case CAF_ATOM_NAME_STR:
        *elem_size = sizeof(char);
        return true;
    case CAF_ATOM_NAMES:
        *elem_size = sizeof(char *);
        return true;
    default:
        return false;
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

CheckpointStatus unserialize_array(const uint8_t *payload, size_t *bytes_read, uint16_t *out_flag,
                                   uint32_t *out_n, void **out_arr)
{
    CheckpointStatus status = validate_array_header(payload, out_flag, out_n);
    if (status != CHECKPOINT_OK) {
        return status;
    }
    *bytes_read = CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(*out_flag) + sizeof(*out_n);

    size_t elem_size = 0;
    if (*out_flag == CAF_ATOM_NAMES) {
        // caller needs to handle this separately since it's an array of strings
        return CHECKPOINT_OK;
    }
    if (!array_flag_elem_size(*out_flag, &elem_size)) {
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
