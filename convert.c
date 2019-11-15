#include <stdint.h>

#define S(i, j) *(src + (i) + (j) * ldsrc)
#define D(i, j) *(dst + (i) + (j) * lddst)

void convert_double_to_float(double *src, uint64_t ldsrc, float *dst,
                             uint64_t lddst, uint64_t m, uint64_t n) {
    uint64_t i, j;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            D(i, j) = (float)S(i, j);
        }
    }
    return;
}

void convert_float_to_double(float *src, uint64_t ldsrc, double *dst,
                             uint64_t lddst, uint64_t m, uint64_t n) {
    uint64_t i, j;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            D(i, j) = (double)S(i, j);
        }
    }
    return;
}
