/**
 * @file Random.c
 * @author Luis Granadillo
 * @brief Random number generator implementation adapted from Numerical Recipes in C for reentrant
 * use and better properties than built-in rand().
 * @version 0.1
 * @date 2026-06-01
 *
 * @copyright Copyright (c) 2026
 *
 */

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7 // TODO: use a C macro, like FLT_EPSILON
#define RNMX (1.0 - EPS)

// ENHANCE: consider implementing a more modern generator, like xoshiro256** or PCG, for better
// performance and statistical properties
// https://en.wikipedia.org/wiki/Xorshift#xoshiro256**

// TODO: create one function for seeding the generator and one for generating random numbers
// both require the state? so it is updated and used reentrantly

/**
 * @brief ran1 from Numerical Recipes in C (Ch. 7.1), adapted to use a pointer for the seed
 *
 * @param idum
 * @return float
 */
float ran1(long *idum)
{
    /* “Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added
    safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
    values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
    successive deviates in a sequence. RNMX should approximate the largest floating value that is
    less than 1. */
    // TODO: change return type to double and update constants accordingly
    int j;
    long k;
    float temp; // TODO: change to double

    // these retain their values between calls to maintain the state of the generator
    // TODO: expose a structure for the generator state for reentrant use
    static long iy = 0;
    static long iv[NTAB];

    // Initialize when iy is zero
    if (*idum <= 0 || !iy) {
        // Be sure to prevent idum = 0.
        if (-(*idum) < 1) {
            *idum = 1;
        } else {
            *idum = -(*idum);
        }

        // Load the shuffle table(after 8 warm-ups).
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0) {
                *idum += IM;
            }
            if (j < NTAB) {
                iv[j] = *idum;
            }
        }
        iy = iv[0];
    }

    // Start here when not initializing.
    k = (*idum) / IQ;

    // Compute idum = (IA * idum) % IM without overflows by Schrage’s method.
    *idum = IA * (*idum - k * IQ) - IR * k;

    if (*idum < 0) {
        *idum += IM;
    }

    // Will be in the range 0..NTAB - 1.
    j = iy / NDIV;

    // Output previously stored value and refill the shuffle table.
    iy = iv[j];
    iv[j] = *idum;
    temp = AM * iy;
    if (temp > RNMX) {
        // Because users don’t expect endpoint values.
        return RNMX;
    } else {
        return temp;
    }
}
