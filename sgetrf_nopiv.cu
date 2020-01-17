#include <stdio.h>
#include <stdint.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "hpl-ai.h"


#define A(i, j) *(A + (i) + (j) * lda)
#define dA(i, j) *(dA + (i) + (j) * ldda)

void sgetrf_nopiv(int m, int n, float *A, int lda) {

    int j;
    int nb = 32;
    int jb = nb;
    int ldda = lda;

    float one = 1.0;
    float none = -1.0;

    float* dA=NULL;

    cublasStatus_t stat;
    cublasHandle_t handle;

    // Use unblock code.
    if( nb > m || nb > n ) {
        sgetrf2_nopiv(m, n, A, lda);
    } else {
        cublasCreate(&handle);
        cublasSetMathMode(handle, CUBLAS_TENSOR_OP_MATH);
        cudaMalloc((void**)&dA, n*ldda*sizeof(float));
        cublasSetMatrix(n, n, sizeof(float), A, lda, dA, ldda);
        int min_mn = m<n ? m : n;
        for(j=0; j<min_mn; j+=nb) {
            if( min_mn - j < nb ) {
                jb = min_mn - j;
            }

            // Factor panel
            if( j!=0 ) {
                cublasGetMatrix(m, jb, sizeof(float), &dA(0, j), ldda, &A(0, j), lda);
            }
            sgetrf2_nopiv(m-j, jb, &A(j, j), lda);

            if( j+jb < n ) {
                cublasSetMatrix(m-j, jb, sizeof(float), &A(j, j), lda, &dA(j, j), ldda);
                cublasStrsm(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, jb, n-j-jb, &one, &dA(j, j), ldda, &dA(j, j+jb), ldda);

                //strsm('L', 'L', 'N', 'U', jb, n-j-jb, 1.0, &A(j, j), lda, &A(j, j+jb), lda);

                if( j+jb < m ) {
                    cublasGemmEx(handle, CUBLAS_OP_N, CUBLAS_OP_N, m-j-jb, n-j-jb, jb,
                            &none, &dA(j+jb, j), CUDA_R_32F, ldda,
                                   &dA(j, j+jb), CUDA_R_32F, ldda,
                            &one, &dA(j+jb, j+jb), CUDA_R_32F, ldda, CUDA_R_32F, CUBLAS_GEMM_DFALT_TENSOR_OP);
                    //cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m-j-jb, n-j-jb, jb,
                    //        &none, &dA(j+jb, j), ldda, &dA(j, j+jb), ldda, &one, &dA(j+jb, j+jb), ldda);
                    //sgemm('N', 'N', m-j-jb, n-j-jb, jb, -1.0, &A(j+jb, j), lda, &A(j, j+jb), lda, 1.0,
                    //     &A(j+jb, j+jb), lda);


                }

            }
        }
        cudaFree(dA);
        cublasDestroy(handle);
    }
    return;
}

void sgetrf2_nopiv(int m, int n, float *A, int lda) {

  int i;

  if (m <= 1 || n == 0) {
    return;
  }

  if (n == 1) {
    for (i = 1; i < m; i++) {
      A(i, 0) /= A(0, 0);
    }
  } else { // Use recursive code

  int n1 = (m > n ? n : m) / 2;
  int n2 = n - n1;

    sgetrf2_nopiv(m, n1, A, lda);

    strsm('L', 'L', 'N', 'U', n1, n2, 1.0, A, lda, &A(0, n1), lda);

    sgemm('N', 'N', m - n1, n2, n1, -1.0, &A(n1, 0), lda, &A(0, n1), lda, 1.0,
          &A(n1, n1), lda);

    sgetrf2_nopiv(m - n1, n2, &A(n1, n1), lda);

  }
  return;
}
