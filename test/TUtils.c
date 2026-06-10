#include "TUtils.h"
#include "CheckpointUtils.h"
#include "State.h"
#include "unity.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// TODO: move to build/test
const char temp_name[] = "temp.log";

void fopen_error(const char *filename, const FILE *file)
{
    if (file == NULL) {
        fprintf(stderr, "Couldn't open file %s: %s\n", filename, strerror(errno));
        TEST_ASSERT_NOT_NULL_MESSAGE(file, "File not opened - check result file");
    }
}

FILE *open_file(const char *filename)
{
    FILE *file = fopen(filename, "rb");
    fopen_error(filename, file);
    return file;
}

void init_temp(FILE **temp_log)
{
    *temp_log = fopen(temp_name, "wb+");
    fopen_error(temp_name, *temp_log);
}

void clean_temp(FILE **temp_log)
{
    int rc = remove(temp_name);
    if (rc) {
        perror("Remove of test log file failed");
    }
    *temp_log = NULL;
}

void close_if_exists(FILE **file)
{
    if (*file) {
        fclose(*file);
        *file = NULL;
    }
}

void assert_file_exists_and_remove(const char *filename)
{
    FILE *f = fopen("file_does_not_exist.txt", "r");
    TEST_ASSERT_NULL_MESSAGE(f, "File should not exist");
    TEST_ASSERT_EQUAL_INT_MESSAGE(ENOENT, errno, "Error should be ENOENT for non-existent file");
    int errno_nexist = errno;

    f = fopen(filename, "rb");

    if (!f) {
        // expect an error from file not existing, so if its a different error, report it
        if (errno != errno_nexist) {
            fprintf(stderr, "Error opening file %s: %s\n", filename, strerror(errno));
            TEST_FAIL_MESSAGE("File should exist but couldn't be opened");
        }
        char error_msg[256];
        sprintf(error_msg, "File %s should exist but couldn't be opened", filename);
        TEST_ASSERT_NOT_NULL_MESSAGE(f, error_msg);
    }

    fclose(f);

    int ret = remove(filename);
    if (ret != 0) {
        fprintf(stderr, "Error deleting file %s: %s\n", filename, strerror(errno));
    }
}

void assert_many_files_exist_and_remove(const char **output_filenames, size_t num_files)
{
    for (size_t i = 0; i < num_files; ++i) {
        assert_file_exists_and_remove(output_filenames[i]);
    }
}

void build_array_header(uint8_t *buffer, uint16_t flag, uint32_t n)
{
    memcpy(buffer, &((const uint16_t){CHECKPOINT_ARRAY_MAGIC}), sizeof(uint16_t));
    memcpy(buffer + sizeof(uint16_t), &flag, sizeof(uint16_t));
    memcpy(buffer + sizeof(uint16_t) + sizeof(uint16_t), &n, sizeof(uint32_t));
}

// TODO: use pack?
static void build_output_format_payload(const OutputFormat *source, OutFormatPayload *dest)
{
    memset(dest, 0, sizeof(*dest));
    dest->type = (uint8_t)source->type;
    dest->is_active = (uint8_t)source->is_active;

    switch (source->type) {
    case OUTPUT_FORMAT_CSV:
        memcpy(dest->data.csv.filename, source->csv.filename, sizeof(dest->data.csv.filename));
        dest->data.csv.schedule.mode = source->csv.schedule.mode;
        dest->data.csv.schedule.interval = source->csv.schedule.interval;
        dest->data.csv.schedule.list_len = source->csv.schedule.list_len;
        dest->data.csv.schedule.list_idx = source->csv.schedule.list_idx;
        dest->data.csv.schedule.next_checkpoint = source->csv.schedule.next_checkpoint;
        dest->data.csv.schedule.frame_num = source->csv.schedule.frame_num;
        dest->data.csv.field_count = source->csv.field_count;
        dest->data.csv.frame_num = source->csv.frame_num;
        break;
    case OUTPUT_FORMAT_XYZ:
        memcpy(dest->data.xyz.prefix, source->xyz.prefix, sizeof(dest->data.xyz.prefix));
        memcpy(dest->data.xyz.suffix, source->xyz.suffix, sizeof(dest->data.xyz.suffix));
        dest->data.xyz.schedule.mode = source->xyz.schedule.mode;
        dest->data.xyz.schedule.interval = source->xyz.schedule.interval;
        dest->data.xyz.schedule.list_len = source->xyz.schedule.list_len;
        dest->data.xyz.schedule.list_idx = source->xyz.schedule.list_idx;
        dest->data.xyz.schedule.next_checkpoint = source->xyz.schedule.next_checkpoint;
        dest->data.xyz.schedule.frame_num = source->xyz.schedule.frame_num;
        dest->data.xyz.frame_num = source->xyz.frame_num;
        dest->data.xyz.stripped = (uint8_t)source->xyz.stripped;
        break;
    case OUTPUT_FORMAT_STEPS_CSV:
        memcpy(dest->data.steps.filename, source->steps.filename,
               sizeof(dest->data.steps.filename));
        dest->data.steps.with_coordination = (uint8_t)source->steps.with_coordination;
        break;
    default:
        TEST_FAIL_MESSAGE("unexpected output format type in round-trip test");
    }
}

static void assert_output_schedule_matches(const OutSchedPayload *expected,
                                           const OutSchedPayload *actual, const char *context)
{
    TEST_ASSERT_EQUAL_INT_MESSAGE(expected->mode, actual->mode, context);
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(expected->interval, actual->interval, context);
    TEST_ASSERT_EQUAL_INT_MESSAGE(expected->list_len, actual->list_len, context);
    TEST_ASSERT_EQUAL_INT_MESSAGE(expected->list_idx, actual->list_idx, context);
    TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(expected->next_checkpoint, actual->next_checkpoint, context);
    TEST_ASSERT_EQUAL_INT_MESSAGE(expected->frame_num, actual->frame_num, context);
}

void assert_output_format_payload_matches(const OutFormatPayload *expected,
                                          const OutFormatPayload *actual, const char *context)
{
    TEST_ASSERT_EQUAL_INT_MESSAGE(expected->type, actual->type, context);
    TEST_ASSERT_EQUAL_INT_MESSAGE(expected->is_active, actual->is_active, context);

    switch (expected->type) {
    case OUTPUT_FORMAT_CSV:
        TEST_ASSERT_EQUAL_STRING_MESSAGE(expected->data.csv.filename, actual->data.csv.filename,
                                         context);
        assert_output_schedule_matches(&expected->data.csv.schedule, &actual->data.csv.schedule,
                                       context);
        TEST_ASSERT_EQUAL_INT_MESSAGE(expected->data.csv.field_count, actual->data.csv.field_count,
                                      context);
        TEST_ASSERT_EQUAL_INT_MESSAGE(expected->data.csv.frame_num, actual->data.csv.frame_num,
                                      context);
        break;
    case OUTPUT_FORMAT_XYZ:
        assert_output_schedule_matches(&expected->data.xyz.schedule, &actual->data.xyz.schedule,
                                       context);
        TEST_ASSERT_EQUAL_STRING_MESSAGE(expected->data.xyz.prefix, actual->data.xyz.prefix,
                                         context);
        TEST_ASSERT_EQUAL_STRING_MESSAGE(expected->data.xyz.suffix, actual->data.xyz.suffix,
                                         context);
        TEST_ASSERT_EQUAL_INT_MESSAGE(expected->data.xyz.frame_num, actual->data.xyz.frame_num,
                                      context);
        TEST_ASSERT_EQUAL_INT_MESSAGE((int)expected->data.xyz.stripped,
                                      (int)actual->data.xyz.stripped, context);
        break;
    case OUTPUT_FORMAT_STEPS_CSV:
        TEST_ASSERT_EQUAL_STRING_MESSAGE(expected->data.steps.filename, actual->data.steps.filename,
                                         context);
        TEST_ASSERT_EQUAL_INT_MESSAGE((int)expected->data.steps.with_coordination,
                                      (int)actual->data.steps.with_coordination, context);
        break;
    default:
        TEST_FAIL_MESSAGE("unexpected output format type in round-trip test");
    }
}

void assert_output_format_matches_runtime_round_trip(const OutputFormat *expected,
                                                     const OutputFormat *actual,
                                                     const char *message)
{
    OutFormatPayload expected_payload = {0};
    build_output_format_payload(expected, &expected_payload);

    OutFormatPayload actual_payload = {0};
    build_output_format_payload(actual, &actual_payload);

    assert_output_format_payload_matches(&expected_payload, &actual_payload, message);
}

void assert_output_format_matches_runtime(const OutputFormat *expected,
                                          const OutFormatPayload *actual, const char *context)
{
    OutFormatPayload expected_payload = {0};
    build_output_format_payload(expected, &expected_payload);
    assert_output_format_payload_matches(&expected_payload, actual, context);
}
