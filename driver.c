#include"hpl-ai.h"
#include<stdio.h>
#include<stdlib.h>

void sgemm( char transa, char transb, uint64_t m, uint64_t n, uint64_t k,
                    float alpha, float* A, uint64_t lda, float* B, uint64_t ldb,
                                float beta, float* C, uint64_t ldc);

void strsm( char side, char uplo, char transa, char diag, uint64_t m,
                    uint64_t n, float alpha, float* A, uint64_t lda, float* B,
                                uint64_t ldb );

int main(int argc, char* argv[]) {

    uint64_t n = 100;
    uint64_t max_it = 50;
    if( argc >= 2 ) {
        n = atoi(argv[1]);
    }
    if( argc >= 3 ) {
        max_it = atoi(argv[2]);
    }
    if( max_it >= n ) {
        max_it = n-1;
    }

    uint64_t lda = (n + 16 - 1) / 16 * 16; // round up to multiple of 16

    double* A = (double*)malloc(lda*n*sizeof(double));
    double* LU = (double*)malloc(lda*n*sizeof(double));
    double* b = (double*)malloc(n*sizeof(double));
    double* x = (double*)malloc(n*sizeof(double));
    float* sA = (float*)malloc(lda*n*sizeof(float));
    float* sb = (float*)malloc(n*sizeof(float));

    matgen(A, lda, n);
    vecgen(b, n);

    convert_double_to_float(A, lda, sA, lda, n, n);
    convert_double_to_float(b, n, sb, n, n, 1);

    sgetrf2_nopiv(n, n, sA, lda);
    
    strsm('L', 'L', 'N', 'U', n, 1, 1.0, sA, lda, sb, n);
    strsm('L', 'U', 'N', 'N', n, 1, 1.0, sA, lda, sb, n);

    convert_float_to_double(sA, lda, LU, lda, n, n);
    convert_float_to_double(sb, n, x, n, n, 1);

    gmres(n, A, lda, x, b, LU, lda, max_it, 1, 1e-15);
    
    double normA = dlange('F', n, n, A, lda);
    double normx = dlange('F', n, 1, x, n);
    dgemv('N', n, n, -1.0, A, lda, x, 1, 1.0, b, 1);
    double error = dlange('F', n, 1, b, n) / normA / normx / n;
    printf("Final backward error |b-A*x| / (|A||x|n) : %e\n",error);

    return 0;
}
