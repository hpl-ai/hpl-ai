#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

int64_t state = 1; // Can be any odd number. Not thread safe.

// Multiplicative Congruential Generator (MCG)
uint32_t mcg_rand() {
  const uint64_t MULTIPLIER = 14647171131086947261U;
  state *= MULTIPLIER;
  return state >> 32;
}

// Generate double floating-point number from uniform(-0.5, 0.5)
double mcg_rand_double() { return ((double)mcg_rand()) / UINT32_MAX - 0.5; }

// Generate a row diagonal dominate square matrix A.
// Column diagonal dominate matrix would also work.
void matgen(double *A, uint64_t lda, uint64_t n) {

  double *diag = (double *)malloc(n * sizeof(double));
  memset(diag, 0, n * sizeof(double));

  for (uint64_t j = 0; j < n; j++) {
    for (uint64_t i = 0; i < n; i++) {
      A[j * lda + i] = mcg_rand_double();
      diag[i] += fabs(A[j * lda + i]);
    }
  }

  // TODO: Use reduction to compute row sum to avoid underflow.
  for (uint64_t i = 0; i < n; i++) {
    A[i * lda + i] = diag[i] - fabs(A[i * lda + i]);
  }

  free(diag);
}

void vecgen(double *v, uint64_t n) {
    for (uint64_t i = 0; i < n; i++) {
        v[i] = mcg_rand_double();
    }
    return;
}
