#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hpl-ai.h"

#define A(i, j) *HPLAI_INDEX2D(A, (i), (j), lda)
#define STATE_INIT 1 // RNG initial state. Can be any odd number.

static const unsigned long int pow32_1 = 4294967295UL;

static unsigned long long int state = STATE_INIT;  // Not thread safe.

// Multiplicative Congruential Generator (MCG)
unsigned long int mcg_rand() {
    const unsigned long long int MULTIPLIER = 14647171131086947261ULL;
    state *= MULTIPLIER;
    return state >> 32; /* use high 32 bits */
}

// Jump ahead function to go through N steps in log(N) time.
void mcg_advance(unsigned int delta) {
    const unsigned long long int MULTIPLIER = 14647171131086947261ULL;
    unsigned long long int accum = MULTIPLIER;
    while(delta != 0) {
        if(delta & 1) {
            delta = delta - 1;
            state *= accum;
        }
        delta = delta / 2;
        accum = accum * accum;
    }
}

// Reset RNG state.
void mcg_reset() {
    state = STATE_INIT;
}

// Generate double floating-point number from uniform(-0.5, 0.5)
double mcg_rand_double() { return ((double)mcg_rand()) / pow32_1 - 0.5; }

// Generate a row diagonally dominant square matrix A.
void matgen(double *A, int lda, int m) {

    int i, j;

    double *diag = (double *)malloc(m * sizeof(double));
    memset(diag, 0, m * sizeof(double));

    for (j = 0; j < m; j++) {
        for (i = 0; i < m; i++) {
            A(i, j) = mcg_rand_double();
            diag[i] += fabs(A(i, j));
        }
    }

    for (i = 0; i < m; i++) {
        A(i, i) = diag[i] - fabs(A(i, i));
    }

    free(diag);
}

void vecgen(double *v, int n) {
    int i;
    for (i = 0; i < n; i++) {
        v[i] = mcg_rand_double();
    }
    return;
}
