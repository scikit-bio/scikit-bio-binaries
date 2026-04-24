/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 5: mat_to_centered correctness under WASM.
 *
 * Uses the same 6x6 unweighted-unifrac distance matrix and the same
 * hard-coded expected centered matrix as the existing native test
 * (src/tests/test_pcoa.cpp:test_center_mat). Tolerance 1e-6 matches
 * the native test's vec_almost_equal.
 *
 * Centering is deterministic, platform-independent arithmetic (no RNG,
 * no BLAS — it's a hand-rolled O(n^2) double-loop in
 * principal_coordinate_analysis.cpp). This test should always pass
 * bit-exactly modulo the accumulated 1e-6 rounding in row_sum.
 */

#include "scikit-bio-binaries/ordination.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

static int failures = 0;
#define CHECK(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, #cond); \
        ++failures; \
    } \
} while (0)

int main() {
    // Input: unweighted-UniFrac distance matrix from test.biom, per
    // src/tests/test_pcoa.cpp:16-26.
    double matrix[] = {
        0.0000000000, 0.2000000000, 0.5714285714, 0.6000000000, 0.5000000000, 0.2000000000,
        0.2000000000, 0.0000000000, 0.4285714286, 0.6666666667, 0.6000000000, 0.3333333333,
        0.5714285714, 0.4285714286, 0.0000000000, 0.7142857143, 0.8571428571, 0.4285714286,
        0.6000000000, 0.6666666667, 0.7142857143, 0.0000000000, 0.3333333333, 0.4000000000,
        0.5000000000, 0.6000000000, 0.8571428571, 0.3333333333, 0.0000000000, 0.6000000000,
        0.2000000000, 0.3333333333, 0.4285714286, 0.4000000000, 0.6000000000, 0.0000000000,
    };
    const unsigned int n = 6;

    // Expected centered matrix from the same native test source.
    const double exp[] = {
         0.05343726,  0.04366213, -0.0329743,  -0.07912698, -0.00495654,  0.01995843,
         0.04366213,  0.073887,    0.04867914, -0.11112434, -0.04973167, -0.00537226,
        -0.0329743,   0.04867914,  0.20714475, -0.07737528, -0.17044974,  0.02497543,
        -0.07912698, -0.11112434, -0.07737528,  0.14830877,  0.11192366,  0.00739418,
        -0.00495654, -0.04973167, -0.17044974,  0.11192366,  0.18664966, -0.07343537,
         0.01995843, -0.00537226,  0.02497543,  0.00739418, -0.07343537,  0.02647959,
    };

    const double tol = 1e-6;

    // fp64 (out-of-place)
    {
        double centered[36];
        skbb_center_distance_matrix_fp64(n, matrix, centered);
        for (int i = 0; i < 36; ++i) {
            CHECK(std::fabs(centered[i] - exp[i]) < tol);
        }
    }

    // fp64 -> fp32 mixed
    {
        float centered_fp32[36];
        skbb_center_distance_matrix_fp64_to_fp32(n, matrix, centered_fp32);
        for (int i = 0; i < 36; ++i) {
            CHECK(std::fabs(static_cast<double>(centered_fp32[i]) - exp[i]) < tol);
        }
    }

    // fp32 out-of-place
    {
        float matrix_fp32[36];
        for (int i = 0; i < 36; ++i) matrix_fp32[i] = static_cast<float>(matrix[i]);
        float centered_fp32[36];
        skbb_center_distance_matrix_fp32(n, matrix_fp32, centered_fp32);
        for (int i = 0; i < 36; ++i) {
            CHECK(std::fabs(static_cast<double>(centered_fp32[i]) - exp[i]) < tol);
        }
    }

    // fp64 in-place
    {
        double buf[36];
        for (int i = 0; i < 36; ++i) buf[i] = matrix[i];
        skbb_center_distance_matrix_fp64(n, buf, buf);
        for (int i = 0; i < 36; ++i) {
            CHECK(std::fabs(buf[i] - exp[i]) < tol);
        }
    }

    std::printf("%s: centering checks, %d failures\n",
                failures == 0 ? "PASS" : "FAIL", failures);
    return failures == 0 ? 0 : 1;
}
