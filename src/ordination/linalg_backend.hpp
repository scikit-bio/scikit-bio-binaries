/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Narrow linear-algebra backend used by PCoA / FSVD.
 *
 * The native build links a BLAS/LAPACKE implementation (OpenBLAS / Netlib).
 * The WebAssembly build ships without one and uses Eigen instead.
 *
 * Backend selection (compile-time):
 *   SKBB_BLAS_BACKEND_EIGEN=1  -> Eigen header-only
 *   otherwise                  -> cblas + LAPACKE
 *
 * The only entry points are the three primitives our PCoA/FSVD pipeline
 * actually needs. We do not attempt a general BLAS facade.
 */

#ifndef SKBB_LINALG_BACKEND_HPP
#define SKBB_LINALG_BACKEND_HPP

#include <stdint.h>

namespace skbb {
namespace linalg {

/*
 * Column-major general matrix multiply: C = A * B.
 * A is (m x k), B is (k x n), C is (m x n).
 * Alpha and beta are fixed at 1.0 and 0.0 respectively (PCoA's only use).
 * All strides equal the leading dimension (no sub-views).
 */
void gemm_nn(uint32_t m, uint32_t n, uint32_t k,
             const double *A, const double *B, double *C);
void gemm_nn(uint32_t m, uint32_t n, uint32_t k,
             const float  *A, const float  *B, float  *C);

/*
 * In-place QR: on input, H is (rows x cols); on output, H holds Q
 * (rows x qcols) where qcols = min(rows, cols). Returns 0 on success.
 *
 * This replaces the {LAPACKE_dgeqrf + LAPACKE_dorgqr} pair used in
 * principal_coordinate_analysis.cpp to build Q from H in-place.
 */
int qr_inplace(uint32_t rows, uint32_t cols, double *H, uint32_t &qcols);
int qr_inplace(uint32_t rows, uint32_t cols, float  *H, uint32_t &qcols);

/*
 * SVD "N, O" variant (matches LAPACKE_dgesvd flags 'N','O'):
 *   - jobu  = 'N' : do not return U
 *   - jobvt = 'O' : overwrite T with V^T
 *
 * Layout contract:
 *   T is a column-major (rows x cols) buffer. On output, V^T occupies
 *   the first k = min(rows,cols) rows of T, i.e. for all (r,c) with
 *   0 <= r < k and 0 <= c < cols, T[r + c*rows] = V^T[r,c]. The
 *   trailing rows in [k, rows) are left undefined by LAPACK when
 *   rows > cols (jobvt='O' only writes the k V^T rows). Consumers MUST
 *   use stride `rows` (the full buffer leading dimension); reading k
 *   packed rows as a (k x cols) contiguous submatrix is invalid and
 *   will observe garbage in the gaps.
 *
 *   S (length >= min(rows,cols)) receives the singular values in
 *   descending order. Returns 0 on success.
 */
int svd_no(uint32_t rows, uint32_t cols, double *T, double *S);
int svd_no(uint32_t rows, uint32_t cols, float  *T, float  *S);

}  // namespace linalg
}  // namespace skbb

#endif
