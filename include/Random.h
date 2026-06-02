#ifndef RANDOM_H
#define RANDOM_H

typedef struct {
    unsigned long long u, v, w;
} RandomState;

void sran(unsigned long long j, RandomState *state);
unsigned long long ran(RandomState *state);
double dran(RandomState *state);

#endif // RANDOM_H
