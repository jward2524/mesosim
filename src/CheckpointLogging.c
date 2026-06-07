#include "CheckpointLogging.h"
#include "CheckpointUtils.h"
#include "State.h"
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void pack_logging_payload(LoggingPayload *payload, const struct LoggingState *ls)
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

static void pack_output_schedule(OutSchedPayload *dest, const OutputSchedule *source)
{
    dest->mode = source->mode;
    dest->interval = source->interval;
    dest->list_len = source->list_len;
    dest->list_idx = source->list_idx;
    dest->next_checkpoint = source->next_checkpoint;
    dest->frame_num = source->frame_num;
}

static CheckpointStatus pack_file_position(FILE *file, fpos_t *pos)
{
    if (!file) {
        // fprintf(stderr, "Checkpoint error: null file pointer was provided for an output format");
        return CHECKPOINT_OK;
    }

    int ret = fgetpos(file, pos);
    if (ret) {
        fprintf(stderr, "Checkpoint error: failed to get file position for CSV format: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }
    return CHECKPOINT_OK;
}

static CheckpointStatus pack_csv_format(OutFormatCsvPayload *dest, const OutputFormat *source)
{
    pack_output_schedule(&dest->schedule, &source->csv.schedule);
    memcpy(dest->filename, source->csv.filename, sizeof(dest->filename));
    dest->field_count = source->csv.field_count;
    dest->frame_num = source->csv.frame_num;
    return pack_file_position(source->csv.file, &dest->file_position);
}

static CheckpointStatus pack_xyz_format(OutFormatXyzPayload *dest, const OutputFormat *source)
{
    pack_output_schedule(&dest->schedule, &source->xyz.schedule);
    memcpy(dest->prefix, source->xyz.prefix, sizeof(dest->prefix));
    memcpy(dest->suffix, source->xyz.suffix, sizeof(dest->suffix));
    dest->frame_num = source->xyz.frame_num;
    dest->stripped = (uint8_t)source->xyz.stripped;
    return CHECKPOINT_OK;
}

static CheckpointStatus pack_steps_format(OutFormatStepsPayload *dest, const OutputFormat *source)
{
    memcpy(dest->filename, source->steps.filename, sizeof(dest->filename));
    dest->with_coordination = (uint8_t)source->steps.with_coordination;
    return pack_file_position(source->steps.file, &dest->file_position);
}

CheckpointStatus pack_output_format_array(OutFormatArrPayload *arr_payload,
                                          const struct LoggingState *ls)
{
    if (arr_payload == NULL || ls == NULL) {
        return CHECKPOINT_ERROR;
    }

    arr_payload->n_out_formats = ls->out_formats_cnt;
    if (arr_payload->n_out_formats <= 0) {
        arr_payload->formats = NULL;
        return CHECKPOINT_OK;
    }

    arr_payload->formats = (OutFormatPayload *)calloc((size_t)arr_payload->n_out_formats,
                                                      sizeof(*arr_payload->formats));
    if (arr_payload->formats == NULL) {
        fprintf(stderr, "Checkpoint error: failed to allocate memory for output formats: %s\n",
                strerror(errno));
        return CHECKPOINT_ERROR;
    }

    for (int i = 0; i < arr_payload->n_out_formats; ++i) {
        const OutputFormat *source_format = &ls->out_formats[i];
        OutFormatPayload *dest_format = &arr_payload->formats[i];

        dest_format->type = (uint8_t)source_format->type;
        dest_format->is_active = (uint8_t)source_format->is_active;

        switch (source_format->type) {
        case OUTPUT_FORMAT_CSV:
            pack_csv_format(&dest_format->data.csv, source_format);
            break;
        case OUTPUT_FORMAT_XYZ:
            pack_xyz_format(&dest_format->data.xyz, source_format);
            break;
        case OUTPUT_FORMAT_STEPS_CSV:
            pack_steps_format(&dest_format->data.steps, source_format);
            break;
        default:
            fprintf(stderr, "Checkpoint error: unsupported output format type %d\n",
                    source_format->type);
            free(arr_payload->formats);
            arr_payload->formats = NULL;
            arr_payload->n_out_formats = 0;
            return CHECKPOINT_ERROR;
        }
    }

    return CHECKPOINT_OK;
}

static CheckpointStatus serialize_fixed_array(uint16_t flag, uint32_t n, const void *arr,
                                              size_t elem_size, uint8_t **p_payload,
                                              uint32_t *p_payload_bytes)
{
    if (p_payload == NULL || p_payload_bytes == NULL) {
        fprintf(stderr, "Checkpoint error: invalid output buffer for array write\n");
        return CHECKPOINT_ERROR;
    }
    if (n > 0 && arr == NULL) {
        fprintf(stderr, "Checkpoint error: array flag %u has non-zero count but NULL data\n",
                (unsigned)flag);
        return CHECKPOINT_ERROR;
    }
    CheckpointStatus status = serialize_array(flag, n, arr, elem_size, p_payload, p_payload_bytes);
    return status;
}

CheckpointStatus serialize_output_format_array(const OutFormatArrPayload *arr_payload,
                                               uint8_t **p_payload, uint32_t *p_payload_bytes)
{
    if (p_payload == NULL || p_payload_bytes == NULL) {
        fprintf(stderr, "Checkpoint error: invalid output buffer for output format array\n");
        return CHECKPOINT_ERROR;
    }

    if (arr_payload == NULL) {
        // nothing to write
        return CHECKPOINT_OK;
    }

    // turns any negative values into zero
    uint32_t n = (arr_payload->n_out_formats > 0) ? (uint32_t)arr_payload->n_out_formats : 0u;
    size_t elem_size = sizeof(OutFormatPayload);

    if (n > 0 && arr_payload->formats == NULL) {
        fprintf(stderr, "Checkpoint error: output format array has non-zero count but NULL data\n");
        return CHECKPOINT_ERROR;
    }

    CheckpointStatus status =
        serialize_fixed_array((uint16_t)CAF_OUTPUT_FORMATS, n, arr_payload->formats, elem_size,
                              p_payload, p_payload_bytes);

    return status;
}

CheckpointStatus unserialize_output_format_array(uint8_t *payload, size_t *total_bytes_read,
                                                 OutFormatArrPayload *arr_payload)
{
    if (payload == NULL || total_bytes_read == NULL || arr_payload == NULL) {
        fprintf(stderr, "Checkpoint error: invalid input for output format array read\n");
        return CHECKPOINT_ERROR;
    }

    uint32_t n = 0;
    uint16_t flag = 0;
    CheckpointStatus status =
        unserialize_array(payload, total_bytes_read, &flag, &n, (void **)&arr_payload->formats);
    if (status != CHECKPOINT_OK) {
        fprintf(stderr, "Checkpoint error: failed to read output format array header\n");
        return status;
    }

    if (flag != CAF_OUTPUT_FORMATS) {
        fprintf(stderr, "Checkpoint error: expected output format array flag %u, got %u\n",
                CAF_OUTPUT_FORMATS, flag);
        free(arr_payload->formats);
        arr_payload->formats = NULL;
        return CHECKPOINT_ERROR;
    }

    arr_payload->n_out_formats = (int)n;
    return CHECKPOINT_OK;
}

void unpack_logging_payload(const LoggingPayload *payload, struct LoggingState *ls)
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

static void unpack_output_schedule(OutputSchedule *dest, const OutSchedPayload *source)
{
    dest->mode = source->mode;
    dest->interval = source->interval;
    dest->list_len = source->list_len;
    dest->list_idx = source->list_idx;
    dest->next_checkpoint = source->next_checkpoint;
    dest->frame_num = source->frame_num;
}

static CheckpointStatus reopen_output_file(const char *filename, const fpos_t *pos, FILE **file)
{
    if (filename[0] == '\0') {
        return CHECKPOINT_OK;
    }

    *file = fopen(filename, "rb+");
    if (!*file) {
        fprintf(stderr, "Failed to open output file %s: %s\n", filename, strerror(errno));
        return CHECKPOINT_ERROR;
    }
    int ret = fsetpos(*file, pos);
    if (ret) {
        fprintf(stderr, "Failed to set file position for output file %s: %s", filename,
                strerror(errno));
        fclose(*file);
        *file = NULL;
        return CHECKPOINT_ERROR;
    }
    return CHECKPOINT_OK;
}

static CheckpointStatus unpack_csv_format(const OutFormatCsvPayload *csv_payload,
                                          OutputFormat *format)
{
    unpack_output_schedule(&format->csv.schedule, &csv_payload->schedule);
    memcpy(format->csv.filename, csv_payload->filename, sizeof(format->csv.filename));
    format->csv.field_count = csv_payload->field_count;
    format->csv.frame_num = csv_payload->frame_num;
    format->csv.field_count = csv_payload->field_count;
    format->csv.frame_num = csv_payload->frame_num;
    return reopen_output_file(format->csv.filename, &csv_payload->file_position, &format->csv.file);
}

static CheckpointStatus unpack_xyz_format(const OutFormatXyzPayload *xyz_payload,
                                          OutputFormat *format)
{
    unpack_output_schedule(&format->xyz.schedule, &xyz_payload->schedule);
    memcpy(format->xyz.prefix, xyz_payload->prefix, sizeof(format->xyz.prefix));
    memcpy(format->xyz.suffix, xyz_payload->suffix, sizeof(format->xyz.suffix));
    format->xyz.frame_num = xyz_payload->frame_num;
    format->xyz.stripped = (bool)xyz_payload->stripped;
    return CHECKPOINT_OK;
}

static CheckpointStatus unpack_steps_format(const OutFormatStepsPayload *steps_payload,
                                            OutputFormat *format)
{
    memcpy(format->steps.filename, steps_payload->filename, sizeof(format->steps.filename));
    format->steps.with_coordination = (bool)steps_payload->with_coordination;
    return reopen_output_file(format->steps.filename, &steps_payload->file_position,
                              &format->steps.file);
}

CheckpointStatus unpack_output_format_array(const OutFormatArrPayload *arr_payload,
                                            struct LoggingState *ls)
{
    if (arr_payload == NULL || ls == NULL) {
        return CHECKPOINT_OK;
    }

    ls->out_formats_cnt = arr_payload->n_out_formats;
    if (ls->out_formats_cnt <= 0) {
        ls->out_formats = NULL;
        return CHECKPOINT_OK;
    }

    ls->out_formats = (OutputFormat *)calloc((size_t)ls->out_formats_cnt, sizeof(*ls->out_formats));
    if (ls->out_formats == NULL) {
        fprintf(stderr,
                "Checkpoint error: failed to allocate memory for output formats in apply: %s\n",
                strerror(errno));
        ls->out_formats_cnt = 0;
        return CHECKPOINT_ERROR;
    }

    for (int i = 0; i < arr_payload->n_out_formats; ++i) {
        const OutFormatPayload *source_format = &arr_payload->formats[i];
        OutputFormat *dest_format = &ls->out_formats[i];

        dest_format->type = (OutputFormatType)source_format->type;
        dest_format->is_active = (bool)source_format->is_active;

        // TODO: re-open files and move fpointers to the correct positions
        CheckpointStatus status;
        switch (dest_format->type) {
        case OUTPUT_FORMAT_CSV:
            status = unpack_csv_format(&source_format->data.csv, dest_format);
            break;
        case OUTPUT_FORMAT_XYZ:
            status = unpack_xyz_format(&source_format->data.xyz, dest_format);
            break;
        case OUTPUT_FORMAT_STEPS_CSV:
            status = unpack_steps_format(&source_format->data.steps, dest_format);
            break;
        default:
            fprintf(stderr, "Checkpoint error: unsupported output format type %d in apply\n",
                    dest_format->type);
            free(ls->out_formats);
            ls->out_formats = NULL;
            ls->out_formats_cnt = 0;
            status = CHECKPOINT_ERROR;
        }

        if (status != CHECKPOINT_OK) {
            fprintf(stderr, "Checkpoint error: failed to unpack output format %d\n", i);
            free(ls->out_formats);
            ls->out_formats = NULL;
            ls->out_formats_cnt = 0;
            return CHECKPOINT_ERROR;
        }
    }
    return CHECKPOINT_OK;
}
