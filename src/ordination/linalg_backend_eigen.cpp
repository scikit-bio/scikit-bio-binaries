/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Header-only Eigen backend used by the WebAssembly build of
 * scikit-bio-binaries, where a system BLAS/LAPACKE is not available.
 *
 * Each routine implements exactly the LAPACK-shaped contract described in
 * linalg_backend.hpp: same dimensions, same in/out buffers, column-major
 * layout. Eigen does the heavy lifting; we just own the storage shuffle.
 */

#include "linalg_backend.hpp"

#include <algorithm>

// Eigen is sensitive to assert triggers in optimized builds; keep them on
// in debug builds only. The WASM build compiles with -O3 -fno-exceptions.
#ifdef NDEBUG
#  define EIGEN_NO_DEBUG
#endif
// We never allocate temporaries inside real-time paths, but our PCoA
// routine is called once per distance matrix; default policy is fine.

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/SVD>

namespace skbb {
namespace linalg {

template <class T>
using ColMat   = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
template <class T>
using CMapConst= Eigen::Map<const ColMat<T>>;
template <class T>
using CMap     = Eigen::Map<ColMat<T>>;

template <class T>
static inline void gemm_nn_T(uint32_t m, uint32_t n, uint32_t k,
                             const T *A, const T *B, T *C) {
    CMapConst<T> a(A, m, k);
    CMapConst<T> b(B, k, n);
    CMap<T>      c(C, m, n);
    c.noalias() = a * b;
}

void gemm_nn(uint32_t m, uint32_t n, uint32_t k,
             const double *A, const double *B, double *C) {
    gemm_nn_T<double>(m, n, k, A, B, C);
}
void gemm_nn(uint32_t m, uint32_t n, uint32_t k,
             const float *A, const float *B, float *C) {
    gemm_nn_T<float>(m, n, k, A, B, C);
}

template <class T>
static inline int qr_inplace_T(uint32_t rows, uint32_t cols, T *H, uint32_t &qcols) {
    qcols = std::min(rows, cols);

    // householderQr() wants an Eigen::Matrix, not a Map; copy in then copy Q back.
    ColMat<T> A(rows, cols);
    {
        CMapConst<T> src(H, rows, cols);
        A = src;
    }

    Eigen::HouseholderQR<ColMat<T>> qr(A);

    // Reconstruct the thin Q as an (rows x qcols) matrix. Eigen's
    // householderQ() object can be multiplied by an identity of the
    // desired size to extract the thin Q explicitly.
    ColMat<T> Q = ColMat<T>::Identity(rows, qcols);
    Q = qr.householderQ() * Q;

    // Write Q back into the caller's H buffer.
    CMap<T> dst(H, rows, qcols);
    dst = Q;
    return 0;
}

int qr_inplace(uint32_t rows, uint32_t cols, double *H, uint32_t &qcols) {
    return qr_inplace_T<double>(rows, cols, H, qcols);
}
int qr_inplace(uint32_t rows, uint32_t cols, float *H, uint32_t &qcols) {
    return qr_inplace_T<float>(rows, cols, H, qcols);
}

template <class T>
static inline int svd_no_T(uint32_t rows, uint32_t cols, T *T_mat, T *S) {
    // Input matrix T is (rows x cols), column-major, in linear buffer.
    ColMat<T> A(rows, cols);
    {
        CMapConst<T> src(T_mat, rows, cols);
        A = src;
    }

    // Match LAPACKE_gesvd(jobu='N', jobvt='O'): compute V^T, no U.
    // ComputeThinV is sufficient because we overwrite the first
    // min(rows,cols) rows of T_mat with V^T.
    Eigen::BDCSVD<ColMat<T>> svd(A, Eigen::ComputeThinV);

    // Singular values, sorted descending (Eigen guarantees this).
    const uint32_t k = std::min(rows, cols);
    for (uint32_t i = 0; i < k; ++i) {
        S[i] = svd.singularValues()(i);
    }

    // V^T overwrite: in LAPACKE with jobvt='O' the first min(m,n) rows of A
    // hold V^T (the rows of V^T are the right singular vectors).
    //
    // Eigen's V matrix is (cols x k); transpose into the first k rows of T.
    ColMat<T> Vt = svd.matrixV().transpose();   // (k x cols)

    // Write Vt into the first k rows of the (rows x cols) T_mat buffer,
    // leaving the trailing rows untouched (the consumer only reads
    // the first k rows, matching LAPACK's behaviour).
    CMap<T> dst(T_mat, rows, cols);
    for (uint32_t r = 0; r < k; ++r) {
        for (uint32_t c = 0; c < cols; ++c) {
            dst(r, c) = Vt(r, c);
        }
    }
    return 0;
}

int svd_no(uint32_t rows, uint32_t cols, double *T, double *S) {
    return svd_no_T<double>(rows, cols, T, S);
}
int svd_no(uint32_t rows, uint32_t cols, float *T, float *S) {
    return svd_no_T<float>(rows, cols, T, S);
}

}  // namespace linalg
}  // namespace skbb
