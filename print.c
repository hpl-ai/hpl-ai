#include <stdio.h>
#include <stdint.h>

#define A(i, j) *(A + (i) + (j) * lda)

void print_matrix_float(float *A, uint64_t lda, uint64_t m, uint64_t n) {
  if (lda < m) {
    return;
  }
  if (m == 1) {
    printf("[ ");
  } else {
    printf("[\n");
  }
  for (uint64_t i = 0; i < m; ++i) {
    for (uint64_t j = 0; j < n; ++j) {
      printf(" %10.6f", A(i, j));
    }
    if (m > 1) {
      printf("\n");
    } else {
      printf(" ");
    }
  }
  printf("];\n");
  return;
}

void print_matrix_double(double *A, uint64_t lda, uint64_t m, uint64_t n) {
  if (lda < m) {
    return;
  }
  if (m == 1) {
    printf("[ ");
  } else {
    printf("[\n");
  }
  for (uint64_t i = 0; i < m; ++i) {
    for (uint64_t j = 0; j < n; ++j) {
      printf(" %14.10e", A(i, j));
    }
    if (m > 1) {
      printf("\n");
    } else {
      printf(" ");
    }
  }
  printf("];\n");
  return;
}
