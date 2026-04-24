/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Native linear-algebra backend: forwards to cblas + LAPACKE. Exact
 * behaviour of the pre-refactor principal_coordinate_analysis.cpp.
 */

#include "linalg_backend.hpp"

#include <algorithm>
#include <cstdlib>
#include <cblas.h>
#include <lapacke.h>

namespace skbb {
namespace linalg {

void gemm_nn(uint32_t m, uint32_t n, uint32_t k,
             const double *A, const double *B, double *C) {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A, m, B, k, 0.0, C, m);
}

void gemm_nn(uint32_t m, uint32_t n, uint32_t k,
             const float  *A, const float  *B, float  *C) {
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0f, A, m, B, k, 0.0f, C, m);
}

int qr_inplace(uint32_t rows, uint32_t cols, double *H, uint32_t &qcols) {
    qcols = std::min(rows, cols);
    double *tau = new double[qcols];
    int rc = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, rows, cols, H, rows, tau);
    if (rc == 0) {
        rc = LAPACKE_dorgqr(LAPACK_COL_MAJOR, rows, qcols, qcols, H, rows, tau);
    }
    delete[] tau;
    return rc;
}

int qr_inplace(uint32_t rows, uint32_t cols, float *H, uint32_t &qcols) {
    qcols = std::min(rows, cols);
    float *tau = new float[qcols];
    int rc = LAPACKE_sgeqrf(LAPACK_COL_MAJOR, rows, cols, H, rows, tau);
    if (rc == 0) {
        rc = LAPACKE_sorgqr(LAPACK_COL_MAJOR, rows, qcols, qcols, H, rows, tau);
    }
    delete[] tau;
    return rc;
}

int svd_no(uint32_t rows, uint32_t cols, double *T, double *S) {
    double *superb = (double *) malloc(sizeof(double) * rows);
    int rc = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'N', 'O', rows, cols,
                            T, rows, S, NULL, rows, NULL, cols, superb);
    free(superb);
    return rc;
}

int svd_no(uint32_t rows, uint32_t cols, float *T, float *S) {
    float *superb = (float *) malloc(sizeof(float) * rows);
    int rc = LAPACKE_sgesvd(LAPACK_COL_MAJOR, 'N', 'O', rows, cols,
                            T, rows, S, NULL, rows, NULL, cols, superb);
    free(superb);
    return rc;
}

}  // namespace linalg
}  // namespace skbb
