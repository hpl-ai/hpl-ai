#include <stdio.h>
#include <stdint.h>
#include "hpl-ai.h"

#define A(i, j) *(A + (i) + (j) * lda)

void sgetrf2_nopiv(uint64_t m, uint64_t n, float *A, uint64_t lda) {

    uint64_t i;

  if (m <= 1 || n == 0) {
    return;
  }

  if (n == 1) {
    for (i = 1; i < m; i++) {
      A(i, 0) /= A(0, 0);
    }
  } else { // Use recursive code

    uint64_t n1 = (m > n ? n : m) / 2;
    uint64_t n2 = n - n1;

    sgetrf2_nopiv(m, n1, A, lda);

    strsm('L', 'L', 'N', 'U', n1, n2, 1.0, A, lda, &A(0, n1), lda);

    sgemm('N', 'N', m - n1, n2, n1, -1.0, &A(n1, 0), lda, &A(0, n1), lda, 1.0,
          &A(n1, n1), lda);

    sgetrf2_nopiv(m - n1, n2, &A(n1, n1), lda);

  }
  return;
}
