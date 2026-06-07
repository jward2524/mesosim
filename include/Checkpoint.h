#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "CheckpointUtils.h"
#include "State.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

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
    uint8_t has_ls;
    // end of header or placeholder for future features
    uint32_t reserved0;
} CheckpointHeader;
#pragma pack(pop)

// ENHANCE: macro for functions exposed only for testing
uint32_t checkpoint_checksum32(const void *data, size_t size);
void set_checkpoint_header_magic(CheckpointHeader *header);
void initialize_checkpoint_header(CheckpointHeader *header);
void finalize_checkpoint_header(CheckpointHeader *header, uint32_t payload_bytes,
                                uint32_t checksum);
bool checkpoint_header_has_valid_magic(const CheckpointHeader *header);
bool checkpoint_header_has_valid_structure(const CheckpointHeader *header);
bool checkpoint_header_is_valid(const CheckpointHeader *header);
bool checkpoint_header_has_valid_checksum(const CheckpointHeader *header, const void *payload,
                                          size_t payload_size);
CheckpointStatus write_checkpoint_header(FILE *file, const CheckpointHeader *header);
CheckpointStatus read_checkpoint_header(FILE *file, CheckpointHeader *header);

CheckpointStatus write_checkpoint_payload(FILE *file, const uint8_t *payload,
                                          const uint32_t payload_bytes);

CheckpointStatus rebuild_rates_and_transitions(struct SimulationState *ss,
                                               struct SimulationEnv *se);

CheckpointStatus verify_payload_size(FILE *file, CheckpointHeader *header);

CheckpointStatus write_checkpoint_file(const char *path, const struct SimulationState *ss,
                                       const struct SimulationEnv *se,
                                       const struct LoggingState *ls);
CheckpointStatus read_checkpoint_file(const char *path, struct SimulationState *ss,
                                      struct SimulationEnv *se, struct LoggingState *ls);

void checkpoint_on_schedule(const struct SimulationState *ss, const struct SimulationEnv *se,
                            struct LoggingState *ls);

#endif // CHECKPOINT_H
