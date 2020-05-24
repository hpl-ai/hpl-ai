
#include"hpl-ai.h"
#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include <math.h>

void
diff(int m, int n, float *sA1, int lda1, float *sA2, int lda2) {
    int i, j;
    float sum=0.0f, x;

    for (i = 0; i < m; ++i) {
        float rowsum = 0.0f;
        for (j = 0; j < n; ++j) {
            x = fabsf( sA1[i+j*lda1] - sA2[i+j*lda2] );
            sum += x;
            rowsum += x;
        }
        printf("%d %e\n", i, rowsum);
    }

    printf("%d %d %e\n", m, n, sum);
}

int
main(int argc, char *argv[]) {
    int n;
    int lda; // round up to multiple of 16

    if (argc <= 1 || sscanf(argv[1], "%d", &n) != 1 || n < 1)
        n=10;

    lda = (n + 16 - 1) / 16 * 16; // round up to multiple of 16
    double* A = (double*)malloc(lda*n*sizeof(double));
    float* sA = (float*)malloc(lda*n*sizeof(float));
    float* sA_gpu = (float*)malloc(lda*n*sizeof(float));

    matgen(A, lda, n, 1313);
    convert_double_to_float(A, lda, sA, lda, n, n);
    convert_double_to_float(A, lda, sA_gpu, lda, n, n);

    sgetrf_nopiv(n, n, sA, lda);
    sgetrf_nopiv_cpu(n, n, sA_gpu, lda);

    diff(n, n, sA, lda, sA_gpu, lda);

    return 0;
}
