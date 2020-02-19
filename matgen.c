#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hpl-ai.h"

#define A(i, j) *HPLAI_INDEX2D(A, (i), (j), lda)

#define LCG_STATE_INIT 1
#define LCG_MUL 5.4210108624275222e-20f
#define LCG_A  6364136223846793005ULL
#define LCG_C  1ULL

#define MCG_STATE_INIT 1 // RNG initial state. Can be any odd number.
#define MCG_A 14647171131086947261ULL
#define MCG_MUL 2.328306436538696e-10; // 2^-32

// Multiplicative Congruential Generator (MCG)
inline unsigned long int mcg_rand(unsigned long long *piseed) {
    *piseed *= MCG_A;
    return *piseed >> 32; /* use high 32 bits */
}

// Jump ahead function to go through N steps in log(N) time.
inline void mcg_advance(unsigned int delta, unsigned long long* piseed) {
    unsigned long long int accum = MCG_A;
    while(delta != 0) {
        if(delta & 1) {
            delta = delta - 1;
            *piseed *= accum;
        }
        delta = delta / 2;
        accum = accum * accum;
    }
}

// Generate double floating-point number from uniform(-0.5, 0.5)
inline double mcg_rand_double(unsigned long long* piseed) {
    return ((double)mcg_rand(piseed)) * MCG_MUL - 0.5;
}

// Generate a row diagonally dominant square matrix A.
void matgen(double *A, int lda, int m, unsigned long long iseed) {

    int i, j;

    double *diag = (double *)malloc(m * sizeof(double));
    memset(diag, 0, m * sizeof(double));

    for (j = 0; j < m; j++) {
        for (i = 0; i < m; i++) {
            A(i, j) = mcg_rand_double(&iseed);
            diag[i] += fabs(A(i, j));
        }
    }

    for (i = 0; i < m; i++) {
        A(i, i) = diag[i] - fabs(A(i, i));
    }

    free(diag);
}

void vecgen(double *v, int n, unsigned long long iseed) {
    int i;
    for (i = 0; i < n; i++) {
        v[i] = mcg_rand_double(&iseed);
    }
    return;
}
