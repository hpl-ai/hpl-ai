#include"hpl-ai.h"
#include<stdio.h>
#include<stdlib.h>
#include<float.h>

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

    double time_convert, time_factor, time_solve, time_gmres, time_total;

    uint64_t lda = (n + 16 - 1) / 16 * 16; // round up to multiple of 16

    double* A = (double*)malloc(lda*n*sizeof(double));
    double* LU = (double*)malloc(lda*n*sizeof(double));
    double* b = (double*)malloc(n*sizeof(double));
    double* x = (double*)malloc(n*sizeof(double));
    float* sA = (float*)malloc(lda*n*sizeof(float));
    float* sb = (float*)malloc(n*sizeof(float));

    matgen(A, lda, n);
    vecgen(b, n);

    // Convert A and b to single.
    time_convert = get_wtime();
    time_total = time_convert;
    convert_double_to_float(A, lda, sA, lda, n, n);
    convert_double_to_float(b, n, sb, n, n, 1);
    time_convert = get_wtime() - time_convert;
    printf("Time spent in conversion to single: %e second.\n", time_convert);

    // LU factorization without pivoting.
    time_factor = get_wtime();
    time_total = time_factor;
    sgetrf2_nopiv(n, n, sA, lda);
    time_factor = get_wtime() - time_factor;
    printf("Time spent in factorization       : %e second.\n", time_factor);
    
    // Forward and backward substitution.
    time_solve = get_wtime();
    strsm('L', 'L', 'N', 'U', n, 1, 1.0, sA, lda, sb, n);
    strsm('L', 'U', 'N', 'N', n, 1, 1.0, sA, lda, sb, n);
    time_solve = get_wtime() - time_solve;
    printf("Time spent in solve               : %e second.\n", time_solve);

    // Convert result back to double.
    time_convert = get_wtime();
    convert_float_to_double(sA, lda, LU, lda, n, n);
    convert_float_to_double(sb, n, x, n, n, 1);
    time_convert = get_wtime() - time_convert;
    printf("Time spent in conversion to double: %e second.\n", time_convert);

    // Using GMRES without restart.
    time_gmres = get_wtime();
    gmres(n, A, lda, x, b, LU, lda, max_it, 1, 1e-15);
    time_gmres = get_wtime() - time_gmres;
    time_total = get_wtime() - time_total;
    printf("Time spent in GMRES               : %e second.\n", time_gmres);
    printf("Total time                        : %e second.\n", time_total);

    double ops = 2.0/3.0 * n * n * n + 3.0/2.0 * n * n;
    printf("Effective operation per sec       : %12f GFLOPs\n", 1e-9 * ops / time_total );

    
    // Chck final backward error.
    double norm_A = dlange('F', n, n, A, lda);
    double norm_x = dlange('F', n, 1, x, n);
    double norm_b = dlange('F', n, 1, b, n);
    dgemv('N', n, n, 1.0, A, lda, x, 1, -1.0, b, 1);
    double error = dlange('F', n, 1, b, n) / (norm_A * norm_x + norm_b) / n;
    printf("Final backward error |A*x-b| / (|A||x|+|b|) / n : %e\n", error);
    printf("Machine epsilon in double precision: %e\n", DBL_EPSILON);
    if( error < DBL_EPSILON ) {
        printf("Backward error PASSED the machine epsilon threshold.\n");
    } else {
        printf("Backward error DID NOT PASS the macghine epsilon threshold.\n");
    }

    return 0;
}
