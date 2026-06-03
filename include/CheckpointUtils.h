#ifndef CHECKPOINT_UTILS_H
#define CHECKPOINT_UTILS_H

#include <stddef.h>
#include <stdint.h>

#define CHECKPOINT_MAGIC_STRING "MESOSIM2"
#define CHECKPOINT_MAGIC_SIZE 8u
#define CHECKPOINT_FORMAT_VERSION 1u
#define CHECKPOINT_ARRAY_MAGIC 0x4D32 // 'M2' in hex
#define CHECKPOINT_ARRAY_MAGIC_SIZE 2u

typedef enum {
    CHECKPOINT_OK = 0,
    CHECKPOINT_UNIMPLEMENTED = -1,
    CHECKPOINT_ERROR = -2
} CheckpointStatus;

typedef enum ArrFlag {
    CAF_ATOMS = 1,
    CAF_SUBSTRATE_COMPOSITION,
    // CAF_ATOMS_PER_NN_LEVEL,
    CAF_NN_ENERGY,
    CAF_IS_SOLUBLE,
    CAF_ATOM_NAMES,
    CAF_ATOM_NAME_STR,
    CAF_OUTPUT_FORMATS,
    CAF_SENTINEL,
} ArrFlag;

CheckpointStatus append_to_payload(const void *data, size_t data_size, uint8_t **p_payload,
                                   uint32_t *p_payload_bytes);
CheckpointStatus serialize_array_magic(uint8_t **p_payload, uint32_t *p_payload_bytes);
CheckpointStatus validate_array_header(const uint8_t *header, uint16_t *out_flag, uint32_t *out_n);
CheckpointStatus serialize_array_header_into_payload(uint16_t flag, uint32_t n,
                                                     uint8_t **p_payload);
CheckpointStatus serialize_array(uint16_t flag, uint32_t n, const void *arr, size_t elem_size,
                                 uint8_t **p_payload, uint32_t *p_payload_bytes);
CheckpointStatus unserialize_array(const uint8_t *payload, size_t *bytes_read, uint16_t *out_flag,
                                   uint32_t *out_n, void **out_arr);

#endif // CHECKPOINT_UTILS_H
