#include "TUtils.h"
#include "State.h"
#include "unity.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    FILE *file = fopen(filename, "r");
    fopen_error(filename, file);
    return file;
}

void init_temp(FILE **temp_log)
{
    *temp_log = fopen(temp_name, "w+");
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
