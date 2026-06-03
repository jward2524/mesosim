#include "CheckpointUtils.h"
#include "TUtils.h"
#include "unity.h"
#include <stdlib.h>
#include <string.h>

void setUp(void)
{
}

void tearDown(void)
{
}

void test_read_array_parses_double_array_payload(void)
{
    // Parse a numeric array and verify the helper reports the bytes it consumed.
    const double expected[] = {1.25, -2.5, 3.75};
    uint8_t payload[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t) +
                    sizeof(expected)] = {0};
    uint16_t out_flag = 0u;
    uint32_t out_n = 0u;
    size_t bytes_read = 0u;
    double *out_arr = NULL;

    build_array_header(payload, CAF_NN_ENERGY, (uint32_t)(sizeof(expected) / sizeof(expected[0])));
    memcpy(payload + CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t), expected,
           sizeof(expected));

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_OK,
        unserialize_array(payload, &bytes_read, &out_flag, &out_n, (void **)&out_arr),
        "valid numeric array payload should parse successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(CAF_NN_ENERGY, out_flag,
                                   "parsed array should report the expected flag");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(3u, out_n, "parsed array should report the expected length");
    TEST_ASSERT_EQUAL_UINT_MESSAGE((unsigned int)(CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) +
                                                  sizeof(uint32_t) + sizeof(expected)),
                                   (unsigned int)bytes_read,
                                   "parsed array should report the consumed byte count");
    TEST_ASSERT_NOT_NULL_MESSAGE(out_arr, "parsed array should allocate output storage");
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY_MESSAGE(expected, out_arr, 3,
                                           "parsed array should round-trip element values");

    free(out_arr);
}

void test_read_array_rejects_corrupted_magic_without_consuming_bytes(void)
{
    // Corrupt the array magic and confirm the helper leaves the output parameters untouched.
    uint8_t payload[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    uint16_t out_flag = 0u;
    uint32_t out_n = 0u;
    size_t bytes_read = 123u;
    void *out_arr = (void *)0x1;

    build_array_header(payload, CAF_SUBSTRATE_COMPOSITION, 2u);
    payload[0] ^= 0xFFu;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_ERROR, unserialize_array(payload, &bytes_read, &out_flag, &out_n, &out_arr),
        "corrupted array magic should be rejected");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(123u, (unsigned int)bytes_read,
                                   "failed parse should not advance the byte count");
    TEST_ASSERT_EQUAL_PTR_MESSAGE((void *)0x1, out_arr,
                                  "failed parse should leave the output pointer untouched");
}

void test_write_array_header_into_payload_writes_expected_header_bytes(void)
{
    // Write a standalone array header into a preallocated buffer and inspect the bytes directly.
    uint8_t buffer[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    uint8_t *buffer_ptr = buffer;
    const uint16_t expected_flag = CAF_ATOM_NAMES;
    const uint32_t expected_n = 5u;

    CheckpointStatus status =
        serialize_array_header_into_payload(expected_flag, expected_n, &buffer_ptr);
    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, status,
                                  "valid output buffer should accept a header write");

    // Build an expected magic value in-place so the test can compare against raw bytes.
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&((const uint16_t){CHECKPOINT_ARRAY_MAGIC}), buffer,
                                     sizeof(uint16_t),
                                     "written header should start with the array magic");

    TEST_ASSERT_EQUAL_UINT16_MESSAGE(expected_flag, *(uint16_t *)(buffer + sizeof(uint16_t)),
                                     "written header should encode the array flag");
    TEST_ASSERT_EQUAL_UINT32_MESSAGE(expected_n,
                                     *(uint32_t *)(buffer + sizeof(uint16_t) + sizeof(uint16_t)),
                                     "written header should encode the element count");
}

void test_write_array_header_into_payload_rejects_null_output_buffer(void)
{
    // A missing output buffer should fail immediately.
    uint8_t *buffer_ptr = NULL;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR,
                                  serialize_array_header_into_payload(CAF_ATOMS, 1u, &buffer_ptr),
                                  "null output storage should be rejected");
}

void test_write_array_serializes_double_array_payload(void)
{
    // Serialize a numeric array and compare the complete byte layout against the expected bytes.
    const double expected_values[] = {1.5, -2.0, 7.25};
    const size_t header_size = CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t);
    const size_t expected_size = header_size + sizeof(expected_values);
    uint8_t expected_payload[header_size + sizeof(expected_values)];
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;

    build_array_header(expected_payload, CAF_NN_ENERGY,
                       (uint32_t)(sizeof(expected_values) / sizeof(expected_values[0])));
    memcpy(expected_payload + header_size, expected_values, sizeof(expected_values));

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_OK,
        serialize_array(CAF_NN_ENERGY,
                        (uint32_t)(sizeof(expected_values) / sizeof(expected_values[0])),
                        expected_values, sizeof(expected_values[0]), &payload, &payload_bytes),
        "valid array data should serialize successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE((unsigned int)expected_size, payload_bytes,
                                   "serialized array should report the full byte count");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(expected_payload, payload, expected_size,
                                     "serialized array should match the expected byte layout");

    free(payload);
}

void test_write_array_rejects_null_payload_destination(void)
{
    // The helper should reject a missing destination buffer before allocating anything.
    const double values[] = {3.0, 4.0};
    uint32_t payload_bytes = 0u;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_ERROR,
        serialize_array(CAF_NN_ENERGY, 2u, values, sizeof(values[0]), NULL, &payload_bytes),
        "null payload destination should be rejected");
}

void test_append_to_payload_appends_bytes_to_existing_payload(void)
{
    // Start with a tiny existing buffer, then append bytes and verify the old content stays intact.
    uint8_t *payload = (uint8_t *)malloc(2u);
    uint32_t payload_bytes = 2u;
    const uint8_t append_bytes[] = {0x30u, 0x40u, 0x50u};
    const uint8_t expected_bytes[] = {0x10u, 0x20u, 0x30u, 0x40u, 0x50u};

    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "payload should allocate");
    payload[0] = 0x10u;
    payload[1] = 0x20u;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_OK,
        append_to_payload(append_bytes, sizeof(append_bytes), &payload, &payload_bytes),
        "append should succeed with valid inputs");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(5u, payload_bytes,
                                   "payload byte count should include the appended bytes");
    TEST_ASSERT_EQUAL_UINT8_ARRAY_MESSAGE(expected_bytes, payload, 5u,
                                          "payload contents should preserve the original bytes");

    free(payload);
}

void test_append_to_payload_rejects_null_input_buffer(void)
{
    // Passing data without a source buffer should be treated as an immediate error.
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        CHECKPOINT_ERROR, append_to_payload((const uint8_t[]){0xAAu}, 1u, &payload, &payload_bytes),
        "null payload input should be rejected");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(NULL, payload,
                                  "failed append should leave the payload pointer unchanged");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(0u, payload_bytes,
                                   "failed append should leave the payload size unchanged");
}

void test_write_array_magic_writes_expected_magic_bytes(void)
{
    // Verify the helper appends just the two-byte array magic and updates the byte count.
    uint8_t *payload = NULL;
    uint32_t payload_bytes = 0u;
    const uint16_t expected_magic = CHECKPOINT_ARRAY_MAGIC;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, serialize_array_magic(&payload, &payload_bytes),
                                  "writing array magic should succeed with valid inputs");
    TEST_ASSERT_NOT_NULL_MESSAGE(payload, "array magic write should allocate payload bytes");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(2u, payload_bytes, "array magic write should append two bytes");
    TEST_ASSERT_EQUAL_MEMORY_MESSAGE(&expected_magic, payload, sizeof(expected_magic),
                                     "array magic bytes should match the expected value");

    free(payload);
}

void test_write_array_magic_rejects_null_length_pointer(void)
{
    // A null byte-count pointer should fail before any allocation or mutation occurs.
    uint8_t *payload = NULL;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR, serialize_array_magic(&payload, NULL),
                                  "null byte-count storage should be rejected");
    TEST_ASSERT_EQUAL_PTR_MESSAGE(NULL, payload,
                                  "failed write should leave the payload pointer unchanged");
}

void test_validate_array_header_accepts_valid_header(void)
{
    // The raw header bytes should decode cleanly when the magic, flag, and length match.
    uint8_t header[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    uint16_t out_flag = 0u;
    uint32_t out_n = 0u;

    build_array_header(header, CAF_NN_ENERGY, 3u);

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_OK, validate_array_header(header, &out_flag, &out_n),
                                  "valid header should parse successfully");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(CAF_NN_ENERGY, out_flag,
                                   "valid header should report the expected flag");
    TEST_ASSERT_EQUAL_UINT_MESSAGE(3u, out_n, "valid header should report the expected length");
}

void test_validate_array_header_rejects_invalid_magic(void)
{
    // Flip the magic bytes to make sure the header validator rejects a malformed record.
    uint8_t header[CHECKPOINT_ARRAY_MAGIC_SIZE + sizeof(uint16_t) + sizeof(uint32_t)] = {0};
    uint16_t out_flag = 0xFFFFu;
    uint32_t out_n = 0xFFFFFFFFu;

    build_array_header(header, CAF_IS_SOLUBLE, 2u);
    header[0] ^= 0xFFu;

    TEST_ASSERT_EQUAL_INT_MESSAGE(CHECKPOINT_ERROR,
                                  validate_array_header(header, &out_flag, &out_n),
                                  "corrupted magic should be rejected");
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_read_array_parses_double_array_payload);
    RUN_TEST(test_read_array_rejects_corrupted_magic_without_consuming_bytes);

    RUN_TEST(test_validate_array_header_accepts_valid_header);
    RUN_TEST(test_validate_array_header_rejects_invalid_magic);

    RUN_TEST(test_write_array_header_into_payload_writes_expected_header_bytes);
    RUN_TEST(test_write_array_header_into_payload_rejects_null_output_buffer);
    RUN_TEST(test_write_array_serializes_double_array_payload);
    RUN_TEST(test_write_array_rejects_null_payload_destination);

    UNITY_END();
    return 0;
}
