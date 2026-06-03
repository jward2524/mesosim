#ifndef TUTIL_H
#define TUTIL_H

#include "CheckpointLogging.h"
#include <stdint.h>
#include <stdio.h>

#define EXPECT_EXIT(exit_code, block)                                                              \
    do {                                                                                           \
        expected_exit_errno = (exit_code);                                                         \
        if (setjmp(test_exit_jmp) == 0) {                                                          \
            jmp_set = 1;                                                                           \
            block                                                                                  \
        } else {                                                                                   \
            TEST_ASSERT_EQUAL_INT_MESSAGE(expected_exit_errno, exit_errno,                         \
                                          "Expected exit errno does not match actual exit errno"); \
        }                                                                                          \
    } while (0)

void fopen_error(const char *filename, const FILE *file);
FILE *open_file(const char *filename);
void init_temp(FILE **temp_log);
void clean_temp(FILE **temp_log);
void close_if_exists(FILE **file);

void build_array_header(uint8_t *buffer, uint16_t flag, uint32_t n);
void assert_output_format_matches_runtime(const OutputFormat *expected,
                                          const OutFormatPayload *actual, const char *context);
void assert_output_format_payload_matches(const OutFormatPayload *expected,
                                          const OutFormatPayload *actual, const char *context);
void assert_output_format_matches_runtime_round_trip(const OutputFormat *expected,
                                                     const OutputFormat *actual,
                                                     const char *message);

#endif // TUTIL_H
