#pragma once
#include <stdint.h>

void matgen(double *A, uint64_t lda, uint64_t n);
void vecgen(double *v, uint64_t n);
double get_wtime( void );
void print_matrix_float(float *A, uint64_t lda, uint64_t m, uint64_t n);
void print_matrix_double(double *A, uint64_t lda, uint64_t m, uint64_t n);
void convert_double_to_float(double *src, uint64_t ldsrc, float *dst,
                             uint64_t lddst, uint64_t m, uint64_t n);
void convert_float_to_double(float *src, uint64_t ldsrc, double *dst,
                             uint64_t lddst, uint64_t m, uint64_t n);
void sgetrf_nopiv(uint64_t m, uint64_t n, float *A, uint64_t lda);
void sgetrf2_nopiv(uint64_t m, uint64_t n, float *A, uint64_t lda);
void gmres(uint64_t n, double* A, uint64_t lda, double* x, double* b, double* LU, uint64_t ldlu, uint64_t restart, uint64_t max_it, double tol );


//BLAS

void sgemm( char transa, char transb, uint64_t m, uint64_t n, uint64_t k,
                            float alpha, float* A, uint64_t lda, float* B, uint64_t ldb,
                                                            float beta, float* C, uint64_t ldc);

void strsm( char side, char uplo, char transa, char diag, uint64_t m,
                            uint64_t n, float alpha, float* A, uint64_t lda, float* B,
                                                            uint64_t ldb );

void dtrsm( char side, char uplo, char transa, char diag, uint64_t m,
                            uint64_t n, double alpha, double* A, uint64_t lda, double* B,
                                                            uint64_t ldb );


void dgemv(char trans, uint64_t m, uint64_t n, double alpha, double* A, uint64_t lda, double* X, uint64_t incx, double beta, double* Y, uint64_t incy);

double dlange(char norm, uint64_t m, uint64_t n, double* A, uint64_t lda);
