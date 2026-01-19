#ifndef TUTIL_H
#define TUTIL_H

#include <stdio.h>

void fopen_error(const char *filename, const FILE *file);
void init_temp(FILE **temp_log);
void clean_temp(FILE **temp_log);
void close_if_exists(FILE **file);

#endif // TUTIL_H
