#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

int64_t state = 1;  // Can be any odd number. Not thread safe.

// Multiplicative Congruential Generator (MCG)
uint32_t mcg_rand() {
    const unsigned long long int MULTIPLIER = 14647171131086947261ULL;
    state *= MULTIPLIER;
    return state >> 32;
}

// Generate double floating-point number from uniform(-0.5, 0.5)
double mcg_rand_double() { return ((double)mcg_rand()) / UINT32_MAX - 0.5; }

// Generate a row diagonally dominant square matrix A.
void matgen(double *A, int lda, int m) {

    int i, j;

    double *diag = (double *)malloc(m * sizeof(double));
    memset(diag, 0, m * sizeof(double));

    for (j = 0; j < m; j++) {
        for (i = 0; i < m; i++) {
            A[j * lda + i] = mcg_rand_double();
            diag[i] += fabs(A[j * lda + i]);
        }
    }

    for (i = 0; i < m; i++) {
        A[i * lda + i] = diag[i] - fabs(A[i * lda + i]);
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
