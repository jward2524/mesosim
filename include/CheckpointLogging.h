#ifndef CHECKPOINT_LOGGING_H
#define CHECKPOINT_LOGGING_H

#include "CheckpointUtils.h"
#include "State.h"
#include <stdint.h>

typedef struct {
    int framenum;
    int verbose;
    unsigned long verbose_interval;
    int increment_precision;
    int stime_precision;
    int overpot_precision;
} LoggingPayload;

void pack_logging_payload(LoggingPayload *payload, const struct LoggingState *ls);
void unpack_logging_payload(const LoggingPayload *payload, struct LoggingState *ls);

#pragma pack(push, 1)
typedef struct {
    int mode;
    double interval;
    int list_len;
    int list_idx;
    double next_checkpoint;
    int frame_num;
} OutSchedPayload;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct {
    OutSchedPayload schedule;
    char filename[256];
    int field_count;
    int frame_num;
} OutFormatCsvPayload;
typedef struct {
    OutSchedPayload schedule;
    char prefix[256];
    char suffix[256];
    int frame_num;
    uint8_t stripped;
} OutFormatXyzPayload;
typedef struct {
    char filename[256];
    uint8_t with_coordination;
} OutFormatStepsPayload;
typedef struct {
    uint8_t type;
    uint8_t is_active;
    union {
        OutFormatCsvPayload csv;
        OutFormatXyzPayload xyz;
        OutFormatStepsPayload steps;
    } data;
} OutFormatPayload;
#pragma pack(pop)

typedef struct {
    OutFormatPayload *formats;
    int n_out_formats;
} OutFormatArrPayload;

/* Output format array helpers */
CheckpointStatus pack_output_format_array(OutFormatArrPayload *arr_payload,
                                          const struct LoggingState *ls);
void unpack_output_format_array(const OutFormatArrPayload *arr_payload, struct LoggingState *ls);
CheckpointStatus serialize_output_format_array(const OutFormatArrPayload *arr_payload,
                                               uint8_t **p_payload, uint32_t *p_payload_bytes);
CheckpointStatus unserialize_output_format_array(uint8_t *payload, size_t *total_bytes_read,
                                                 OutFormatArrPayload *arr_payload);

#endif // CHECKPOINT_LOGGING_H
