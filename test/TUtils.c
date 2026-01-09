#include "TUtils.h"
#include "unity.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>

void fopen_error(const char *filename, FILE *file)
{
    if (file == NULL) {
        fprintf(stderr, "Couldn't open file %s: %s\n", filename, strerror(errno));
        TEST_ASSERT_NOT_NULL_MESSAGE(file, "File not opened - check result file");
    }
}
