#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#define A(i, j) *(A + (i) + (j) * lda)
#define B(i, j) *(B + (i) + (j) * ldb)
#define C(i, j) *(C + (i) + (j) * ldc)

void sgemm(char transa, char transb, uint64_t m, uint64_t n, uint64_t k,
           float alpha, float *A, uint64_t lda, float *B, uint64_t ldb,
           float beta, float *C, uint64_t ldc) {
    uint64_t i, j, l;

    // Only supprt transa=='N', trabsb=='N'
    if (transa != 'N' || transb != 'N') {
        printf("Not supported in SGEMM.\n");
        return;
    }

    if (m == 0 || n == 0) {
        return;
    }

    if ((alpha == 0.0 || k == 0) && beta == 1.0) {
        return;
    }

    if (alpha == 0.0) {
        if (beta == 0.0) {
            for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                    C(i, j) = 0.0;
                }
            }
        } else {
            for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                    C(i, j) = beta * C(i, j);
                }
            }
        }
    }

    for (j = 0; j < n; j++) {
        if (beta == 0.0) {
            for (i = 0; i < m; i++) {
                C(i, j) = 0.0;
            }
        } else {
            for (i = 0; i < m; i++) {
                C(i, j) = beta * C(i, j);
            }
        }
        for (l = 0; l < k; l++) {
            float temp = alpha * B(l, j);
            for (i = 0; i < m; i++) {
                C(i, j) += temp * A(i, l);
            }
        }
    }
    return;
}
void strsm(char side, char uplo, char transa, char diag, uint64_t m, uint64_t n,
           float alpha, float *A, uint64_t lda, float *B, uint64_t ldb) {

    uint64_t i, j, k;

    // Only support side=='L', transa=='N', alpha==1.0.
    if (side != 'L' || transa != 'N' || alpha != 1.0) {
        printf("Not supported in STRSM.\n");
        return;
    }

    if (m == 0 || n == 0) {
        return;
    }

    int nounit = diag == 'N';

    if (uplo == 'U') {
        for (j = 0; j < n; j++) {
            for (k = m - 1; k < m; k--) {  // stop when k=0-1=UINT64_MAX
                if (nounit) {
                    B(k, j) = B(k, j) / A(k, k);
                }
                for (i = 0; i < k; i++) {
                    B(i, j) = B(i, j) - B(k, j) * A(i, k);
                }
            }
        }
    } else {
        for (j = 0; j < n; j++) {
            for (k = 0; k < m; k++) {
                if (nounit) {
                    B(k, j) = B(k, j) / A(k, k);
                }
                for (i = k + 1; i < m; i++) {
                    B(i, j) = B(i, j) - B(k, j) * A(i, k);
                }
            }
        }
    }
    return;
}
void dtrsm(char side, char uplo, char transa, char diag, uint64_t m, uint64_t n,
           double alpha, double *A, uint64_t lda, double *B, uint64_t ldb) {
    uint64_t i, j, k;

    // Only support side=='L', transa=='N', alpha==1.0.
    if (side != 'L' || transa != 'N' || alpha != 1.0) {
        printf("Not supported in DTRSM.\n");
        return;
    }

    if (m == 0 || n == 0) {
        return;
    }

    int nounit = diag == 'N';

    if (uplo == 'U') {
        for (j = 0; j < n; j++) {
            for (k = m - 1; k < m; k--) {  // stop when k=0-1=UINT64_MAX
                if (nounit) {
                    B(k, j) = B(k, j) / A(k, k);
                }
                for (i = 0; i < k; i++) {
                    B(i, j) = B(i, j) - B(k, j) * A(i, k);
                }
            }
        }
    } else {
        for (j = 0; j < n; j++) {
            for (k = 0; k < m; k++) {
                if (nounit) {
                    B(k, j) = B(k, j) / A(k, k);
                }
                for (i = k + 1; i < m; i++) {
                    B(i, j) = B(i, j) - B(k, j) * A(i, k);
                }
            }
        }
    }
    return;
}

double dlange(char norm, uint64_t m, uint64_t n, double *A, uint64_t lda) {
    uint64_t i, j;

    // Frobenius norm
    if (norm == 'F') {
        double sum = 0.0;
        for (j = 0; j < n; ++j) {
            for (i = 0; i < m; ++i) {
                sum += A(i, j) * A(i, j);
            }
        }
        return sqrt(sum);
        // Infinity norm
    } else if (norm == 'I') {
        double *work = (double *)malloc(m * sizeof(double));
        memset(work, 0, m * sizeof(double));
        double max = 0.0;
        for (j = 0; j < n; ++j) {
            for (i = 0; i < m; ++i) {
                work[i] += fabs(A(i, j));
            }
        }
        for (i = 0; i < m; ++i) {
            if (max < work[i]) {
                max = work[i];
            }
        }
        free(work);
        return max;
    }
    return 0;
}

void dgemv(char trans, uint64_t m, uint64_t n, double alpha, double *A,
           uint64_t lda, double *X, uint64_t incx, double beta, double *Y,
           uint64_t incy) {
    uint64_t i, j;
    if (trans != 'N' || incx != 1 || incy != 1) {
        return;
    }

    if (beta != 1.0) {
        if (beta == 0.0) {
            for (i = 0; i < m; ++i) {
                Y[i] = 0;
            }
        } else {
            for (i = 0; i < m; ++i) {
                Y[i] = beta * Y[i];
            }
        }
    }

    if (alpha == 0.0) {
        return;
    }

    for (j = 0; j < n; ++j) {
        double temp = alpha * X[j];
        for (i = 0; i < m; ++i) {
            Y[i] += temp * A(i, j);
        }
    }
    return;
}
