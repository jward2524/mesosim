#include "CheckpointLogging.h"
#include "TUtils.h"
#include "unity.h"
#include <stdlib.h>
#include <string.h>

// TODO: set up ss and se in setUp so it cleans on test fail
void setUp(void)
{
}

void tearDown(void)
{
}

void test_fill_logging_payload_copies_selected_logging_scalars(void)
{
    // Copy a representative logging state into the compact payload and verify every scalar field.
    struct LoggingState ls = {0};
    LoggingPayload payload = {0};

    ls.sim_log = NULL;
    ls.framenum = 17;
    ls.verbose = 1;
    ls.verbose_interval = 250ul;
    ls.increment_precision = 3;
    ls.stime_precision = 4;
    ls.overpot_precision = 5;

    pack_logging_payload(&payload, &ls);

    // Mutate the source after copying so the payload is the only thing under test.
    ls.framenum = 99;

    TEST_ASSERT_EQUAL_INT_MESSAGE(17, payload.framenum, "frame number should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, payload.verbose, "verbose flag should copy into payload");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(250u, (unsigned int)payload.verbose_interval,
                                   "verbose interval should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, payload.increment_precision,
                                  "increment precision should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, payload.stime_precision,
                                  "stime precision should copy into payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(5, payload.overpot_precision,
                                  "overpotential precision should copy into payload");
    TEST_ASSERT_NOT_EQUAL_INT_MESSAGE(ls.framenum, payload.framenum,
                                      "payload should keep the copied frame number");
}

void test_fill_logging_payload_copies_logging_edge_values(void)
{
    // Use a different mix of values to make sure the helper copies each scalar independently.
    struct LoggingState ls = {0};
    LoggingPayload payload = {0};

    ls.framenum = -8;
    ls.verbose = 0;
    ls.verbose_interval = 1000000ul;
    ls.increment_precision = -2;
    ls.stime_precision = 0;
    ls.overpot_precision = 12;

    pack_logging_payload(&payload, &ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(-8, payload.framenum,
                                  "frame number should copy even when negative");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, payload.verbose, "verbose flag should copy exactly");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(1000000u, (unsigned int)payload.verbose_interval,
                                   "verbose interval should copy exactly");
    TEST_ASSERT_EQUAL_INT_MESSAGE(-2, payload.increment_precision,
                                  "increment precision should copy exactly");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, payload.stime_precision,
                                  "stime precision should copy exactly");
    TEST_ASSERT_EQUAL_INT_MESSAGE(12, payload.overpot_precision,
                                  "overpotential precision should copy exactly");
}

void test_apply_logging_payload_restores_selected_logging_scalars(void)
{
    LoggingPayload payload = {0};
    struct LoggingState ls = {0};

    payload.framenum = 17;
    payload.verbose = 1;
    payload.verbose_interval = 250u;
    payload.increment_precision = 3;
    payload.stime_precision = 4;
    payload.overpot_precision = 5;

    unpack_logging_payload(&payload, &ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(17, ls.framenum, "frame number should restore from payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, ls.verbose, "verbose flag should restore from payload");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(250u, (unsigned int)ls.verbose_interval,
                                   "verbose interval should restore from payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(3, ls.increment_precision,
                                  "increment precision should restore from payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(4, ls.stime_precision,
                                  "stime precision should restore from payload");
    TEST_ASSERT_EQUAL_INT_MESSAGE(5, ls.overpot_precision,
                                  "overpotential precision should restore from payload");
    TEST_ASSERT_NULL_MESSAGE(ls.sim_log, "sim_log should not be modified by the helper");
}

void test_fill_output_format_array_payload_serializes_formats(void)
{
    struct LoggingState ls = {0};
    OutFormatArrPayload arr = {0};

    const int n = 4;
    ls.out_formats_cnt = n;
    ls.out_formats = (OutputFormat *)malloc((size_t)n * sizeof(OutputFormat));
    TEST_ASSERT_NOT_NULL_MESSAGE(ls.out_formats, "out_formats should allocate");

    /* CSV format entry */
    ls.out_formats[0].type = OUTPUT_FORMAT_CSV;
    ls.out_formats[0].is_active = true;
    strncpy(ls.out_formats[0].csv.filename, "out.csv", sizeof(ls.out_formats[0].csv.filename));
    ls.out_formats[0].csv.field_count = 2;
    ls.out_formats[0].csv.frame_num = 42;
    ls.out_formats[0].csv.schedule.mode = OUTPUT_SCHEDULE_INTERVAL_ITERATION;
    ls.out_formats[0].csv.schedule.interval = 1.5;
    ls.out_formats[0].csv.schedule.frame_num = 7;

    /* Second CSV format entry */
    ls.out_formats[1].type = OUTPUT_FORMAT_CSV;
    ls.out_formats[1].is_active = false;
    strncpy(ls.out_formats[1].csv.filename, "other.csv", sizeof(ls.out_formats[1].csv.filename));
    ls.out_formats[1].csv.field_count = 4;
    ls.out_formats[1].csv.frame_num = 99;
    ls.out_formats[1].csv.schedule.mode = OUTPUT_SCHEDULE_LIST_TIME;
    ls.out_formats[1].csv.schedule.interval = 0.25;
    ls.out_formats[1].csv.schedule.frame_num = 11;

    /* STEPS format entry */
    ls.out_formats[2].type = OUTPUT_FORMAT_STEPS_CSV;
    ls.out_formats[2].is_active = true;
    strncpy(ls.out_formats[2].steps.filename, "steps.csv",
            sizeof(ls.out_formats[2].steps.filename));
    ls.out_formats[2].steps.with_coordination = true;

    /* XYZ format entry */
    ls.out_formats[3].type = OUTPUT_FORMAT_XYZ;
    ls.out_formats[3].is_active = true;
    strncpy(ls.out_formats[3].xyz.prefix, "prefix", sizeof(ls.out_formats[3].xyz.prefix));
    strncpy(ls.out_formats[3].xyz.suffix, "suffix", sizeof(ls.out_formats[3].xyz.suffix));
    ls.out_formats[3].xyz.frame_num = 13;
    ls.out_formats[3].xyz.stripped = true;
    ls.out_formats[3].xyz.schedule.mode = OUTPUT_SCHEDULE_INTERVAL_TIME;
    ls.out_formats[3].xyz.schedule.interval = 2.5;
    ls.out_formats[3].xyz.schedule.frame_num = 17;

    CheckpointStatus status = pack_output_format_array(&arr, &ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "fill should succeed");
    TEST_ASSERT_EQUAL_INT_MESSAGE(n, arr.n_out_formats, "n_out_formats should copy");

    for (int i = 0; i < n; ++i) {
        char message[64];
        snprintf(message, sizeof(message), "filled output format %d should match", i);
        assert_output_format_matches_runtime(&ls.out_formats[i], &arr.formats[i], message);
    }

    free(ls.out_formats);
    if (arr.formats) {
        free(arr.formats);
    }
}

void test_write_output_format_array_serializes_formats(void)
{
    struct LoggingState ls = {0};
    OutFormatArrPayload arr = {0};
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;

    /* Create three output formats: CSV, STEPS, XYZ to exercise variants */
    ls.out_formats_cnt = 3;
    ls.out_formats = (OutputFormat *)malloc((size_t)ls.out_formats_cnt * sizeof(OutputFormat));
    TEST_ASSERT_NOT_NULL_MESSAGE(ls.out_formats, "out_formats should allocate");

    /* CSV entry */
    ls.out_formats[0].type = OUTPUT_FORMAT_CSV;
    ls.out_formats[0].is_active = true;
    strncpy(ls.out_formats[0].csv.filename, "out.csv", sizeof(ls.out_formats[0].csv.filename));
    ls.out_formats[0].csv.field_count = 2;

    /* STEPS entry */
    ls.out_formats[1].type = OUTPUT_FORMAT_STEPS_CSV;
    ls.out_formats[1].is_active = false;
    strncpy(ls.out_formats[1].steps.filename, "steps_out.csv",
            sizeof(ls.out_formats[1].steps.filename));
    ls.out_formats[1].steps.with_coordination = true;

    /* XYZ entry */
    ls.out_formats[2].type = OUTPUT_FORMAT_XYZ;
    ls.out_formats[2].is_active = true;
    strncpy(ls.out_formats[2].xyz.prefix, "pfx", sizeof(ls.out_formats[2].xyz.prefix));
    strncpy(ls.out_formats[2].xyz.suffix, "sfx", sizeof(ls.out_formats[2].xyz.suffix));
    ls.out_formats[2].xyz.frame_num = 5;
    ls.out_formats[2].xyz.stripped = false;

    CheckpointStatus status = pack_output_format_array(&arr, &ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "fill should succeed");

    status = serialize_output_format_array(&arr, &payload, &payload_bytes);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status,
                                  "serialize_output_format_array should succeed");
    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "write should allocate payload bytes");

    /* At minimum we expect the array magic to be present at the start of the payload. */
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&((const uint16_t){CHECKPOINT_ARRAY_MAGIC}), payload,
                                     sizeof(uint16_t), "payload should start with array magic");

    /* Verify the array header reports the correct number of elements */
    const size_t header_n_offset = sizeof(uint16_t) + sizeof(uint16_t);
    uint32_t reported_n = 0u;
    TEST_ASSERT_TRUE_MESSAGE(payload_bytes >= header_n_offset + sizeof(uint32_t),
                             "payload should contain complete array header");
    memcpy(&reported_n, payload + header_n_offset, sizeof(reported_n));
    TEST_ASSERT_EQUAL_UINT32_MESSAGE((uint32_t)ls.out_formats_cnt, reported_n,
                                     "array header should report correct element count");

    free(payload);
    free(ls.out_formats);
    if (arr.formats) {
        free(arr.formats);
    }
}

void test_read_output_format_array_parses_written_payload(void)
{
    struct LoggingState ls = {0};
    OutFormatArrPayload arr = {0};
    OutFormatArrPayload parsed = {0};
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;
    size_t bytes_read = 0u;

    ls.out_formats_cnt = 2;
    ls.out_formats = (OutputFormat *)malloc((size_t)ls.out_formats_cnt * sizeof(OutputFormat));
    TEST_ASSERT_NOT_NULL_MESSAGE(ls.out_formats, "out_formats should allocate");

    ls.out_formats[0].type = OUTPUT_FORMAT_CSV;
    ls.out_formats[0].is_active = true;
    strncpy(ls.out_formats[0].csv.filename, "read1.csv", sizeof(ls.out_formats[0].csv.filename));
    ls.out_formats[0].csv.field_count = 3;

    ls.out_formats[1].type = OUTPUT_FORMAT_XYZ;
    ls.out_formats[1].is_active = false;
    strncpy(ls.out_formats[1].xyz.prefix, "rpx", sizeof(ls.out_formats[1].xyz.prefix));
    strncpy(ls.out_formats[1].xyz.suffix, "rsx", sizeof(ls.out_formats[1].xyz.suffix));

    CheckpointStatus status = pack_output_format_array(&arr, &ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "fill should succeed");

    status = serialize_output_format_array(&arr, &payload, &payload_bytes);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "write should succeed");
    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "payload should be allocated");

    status = unserialize_output_format_array(payload, &bytes_read, &parsed);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "read should succeed");
    TEST_ASSERT_EQUAL_UINT_MESSAGE((unsigned)arr.n_out_formats, (unsigned)parsed.n_out_formats,
                                   "parsed count should match");

    if (parsed.formats) {
        free(parsed.formats);
    }
    free(payload);
    free(ls.out_formats);
    if (arr.formats) {
        free(arr.formats);
    }
}

void test_read_output_format_array_round_trip(void)
{
    struct LoggingState ls = {0};
    OutFormatArrPayload arr = {0};
    OutFormatArrPayload parsed = {0};
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;
    size_t bytes_read = 0u;

    ls.out_formats_cnt = 2;
    ls.out_formats = (OutputFormat *)malloc((size_t)ls.out_formats_cnt * sizeof(OutputFormat));
    TEST_ASSERT_NOT_NULL_MESSAGE(ls.out_formats, "out_formats should allocate");

    ls.out_formats[0].type = OUTPUT_FORMAT_STEPS_CSV;
    ls.out_formats[0].is_active = true;
    strncpy(ls.out_formats[0].steps.filename, "round_steps.csv",
            sizeof(ls.out_formats[0].steps.filename));
    ls.out_formats[0].steps.with_coordination = false;

    ls.out_formats[1].type = OUTPUT_FORMAT_CSV;
    ls.out_formats[1].is_active = true;
    strncpy(ls.out_formats[1].csv.filename, "round.csv", sizeof(ls.out_formats[1].csv.filename));
    ls.out_formats[1].csv.field_count = 1;

    CheckpointStatus status = pack_output_format_array(&arr, &ls);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "fill should succeed");

    status = serialize_output_format_array(&arr, &payload, &payload_bytes);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "write should succeed");
    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "payload should be allocated");

    status = unserialize_output_format_array(payload, &bytes_read, &parsed);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status, "read should succeed");
    TEST_ASSERT_EQUAL_UINT_MESSAGE((unsigned)arr.n_out_formats, (unsigned)parsed.n_out_formats,
                                   "parsed count should match");

    for (int i = 0; i < arr.n_out_formats; ++i) {
        char message[64];
        snprintf(message, sizeof(message), "output format %d should round-trip", i);
        assert_output_format_payload_matches(&arr.formats[i], &parsed.formats[i], message);
    }

    if (parsed.formats) {
        free(parsed.formats);
    }
    free(payload);
    free(ls.out_formats);
    if (arr.formats) {
        free(arr.formats);
    }
}

void test_apply_output_format_array_payload_restores_formats(void)
{
    struct LoggingState ls = {0};
    OutFormatArrPayload arr = {0};

    arr.n_out_formats = 3;
    arr.formats = (OutFormatPayload *)malloc((size_t)arr.n_out_formats * sizeof(OutFormatPayload));
    TEST_ASSERT_NOT_NULL_MESSAGE(arr.formats, "formats array should allocate");

    arr.formats[0].type = OUTPUT_FORMAT_CSV;
    arr.formats[0].is_active = true;
    strncpy(arr.formats[0].data.csv.filename, "applied.csv",
            sizeof(arr.formats[0].data.csv.filename));
    arr.formats[0].data.csv.field_count = 5;

    arr.formats[1].type = OUTPUT_FORMAT_XYZ;
    arr.formats[1].is_active = false;
    strncpy(arr.formats[1].data.xyz.prefix, "apx", sizeof(arr.formats[1].data.xyz.prefix));
    strncpy(arr.formats[1].data.xyz.suffix, "asx", sizeof(arr.formats[1].data.xyz.suffix));
    arr.formats[1].data.xyz.frame_num = 88;
    arr.formats[1].data.xyz.stripped = true;

    arr.formats[0].type = OUTPUT_FORMAT_STEPS_CSV;
    arr.formats[0].is_active = true;
    strncpy(arr.formats[0].data.steps.filename, "applied_steps.csv",
            sizeof(arr.formats[0].data.steps.filename));
    arr.formats[0].data.steps.with_coordination = 0;

    unpack_output_format_array(&arr, &ls);

    TEST_ASSERT_EQUAL_INT_MESSAGE(2, ls.out_formats_cnt,
                                  "output format count should restore from payload");
    assert_output_format_matches_runtime(&ls.out_formats[0], &arr.formats[0],
                                         "first output format should restore");
    assert_output_format_matches_runtime(&ls.out_formats[1], &arr.formats[1],
                                         "second output format should restore");

    free(arr.formats);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_fill_logging_payload_copies_selected_logging_scalars);
    RUN_TEST(test_fill_logging_payload_copies_logging_edge_values);
    RUN_TEST(test_apply_logging_payload_restores_selected_logging_scalars);
    RUN_TEST(test_fill_output_format_array_payload_serializes_formats);
    RUN_TEST(test_write_output_format_array_serializes_formats);
    RUN_TEST(test_read_output_format_array_parses_written_payload);
    RUN_TEST(test_read_output_format_array_round_trip);

    UNITY_END();
    return 0;
}
