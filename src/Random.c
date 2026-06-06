/**
 * @file Random.c
 * @author Luis Granadillo
 * @brief Random number generator implementation adapted from 'Numerical Recipes in C++' (3rd
 * edition) for reentrant use and better properties than the built-in rand().
 * @version 0.1
 * @date 2026-06-01
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "Random.h"
#include <limits.h>

/**
 * @brief Initialize the random state with a seed value, then run the generator a few times to mix
 * the state.
 *
 * @param seed seed value to initialize the generator
 * @param state pointer to the random state struct to initialize
 */
void sran(unsigned long long seed, RandomState *state)
{
    state->v = 4101842887655102017ULL;
    state->w = 1;

    if (seed == 0) {
        // seed of 0 is not allowed for this generator, so replace it with a large nonzero value
        // in this case, all 1 bits
        seed = ULLONG_MAX;
    }

    // avalanche step from MumurHash3, based on suggestion from
    // www.numerical.recipes/forumarchive/index.php/t-2267.html
    seed ^= seed >> 33;
    seed *= 0xff51afd7ed558ccdULL;
    seed ^= seed >> 33;
    seed *= 0xc4ceb9fe1a85ec53ULL;
    seed ^= seed >> 33;

    state->u = seed ^ state->v;
    ran(state);
    state->v = state->u;
    ran(state);
    state->w = state->v;
    ran(state);
}

/**
 * @brief Combined random number generator, produces a random positive integer ranging from 0 to
 * ULONG_LONG_MAX (usually 2^64-1). Implementation of Ran (Ch. 7.1) from Numerical Recipes in C++
 * (3rd edition), adapted for reentrant use. Described by [A1_l(C3) + A3_r] ^ B1
 *
 * @param state pointer to the random state struct
 * @return random integer in the range [0, ULONG_LONG_MAX]
 */
unsigned long long ran(RandomState *state)
{
    // u: LCG moduluo 2^64 with A-C set C3 (Ch. 7.1.2 C)
    state->u = state->u * 2862933555777941757LL + 7046029254386353087LL;

    // v: 64-bit xorshift with A3 triples (Ch. 7.1.2 A)
    state->v ^= state->v >> 17;
    state->v ^= state->v << 31;
    state->v ^= state->v >> 8;

    // w: multiply-with-carry with multiplier B1 (Ch. 7.1.2 B) and 32-bit base
    state->w = 4294957665U * (state->w & 0xffffffff) + (state->w >> 32);

    // composed: 64-bit xorshift with A1 triples (Ch. 7.1.2 A) on u
    // x = A1(C3)
    unsigned long long x = state->u ^ (state->u << 21);
    x ^= x >> 35;
    x ^= x << 4;

    // combine the outputs of the three generators with + and ^ (Ch. 7.1.3)
    return (x + state->v) ^ state->w;
}

/**
 * @brief Return a random floating point number in the range [0, 1) using the ran() generator
 *
 * @param state pointer to the random state struct
 * @return random floating-point number in the range [0, 1)
 */
double dran(RandomState *state)
{
    // divide by 2^64 + 1 to get a double in the range [0, 1)
    // division is expensive relative to multiplication, so multiply by the pre-computed reciprocal
    // of 2^64 (written with 18 significant digits, even though the limit of double precision is ~16
    // digits)
    return 5.42101086242752217e-20 * (double)ran(state);
}
