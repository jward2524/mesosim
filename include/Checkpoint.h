#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "State.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

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

// pack/pop makes the struct have 1-byte alignment instead of default (4?)
#pragma pack(push, 1)
typedef struct {
    char magic[CHECKPOINT_MAGIC_SIZE];
    uint32_t format_version;
    uint32_t header_bytes;
    uint32_t payload_bytes;
    uint32_t checksum;
    uint8_t has_ss;
    uint8_t has_se;
    uint8_t has_atoms;
    // end of header or placeholder for future features
    uint32_t reserved0;
} CheckpointHeader;
#pragma pack(pop)

CheckpointStatus checkpoint_save(const char *path, const struct SimulationState *ss,
                                 const struct SimulationEnv *se);
CheckpointStatus checkpoint_load(const char *path, struct SimulationState *ss,
                                 struct SimulationEnv *se);

uint32_t checkpoint_checksum32(const void *data, size_t size);
void checkpoint_header_set_magic(CheckpointHeader *header);
void checkpoint_header_init(CheckpointHeader *header);
void checkpoint_header_finalize(CheckpointHeader *header, uint32_t payload_bytes,
                                uint32_t checksum);
bool checkpoint_header_has_valid_magic(const CheckpointHeader *header);
bool checkpoint_header_has_valid_structure(const CheckpointHeader *header);
bool checkpoint_header_has_valid_checksum(const CheckpointHeader *header, const void *payload,
                                          size_t payload_size);
CheckpointStatus checkpoint_header_write(FILE *file, const CheckpointHeader *header);
CheckpointStatus checkpoint_header_read(FILE *file, CheckpointHeader *header);

#endif // CHECKPOINT_H
