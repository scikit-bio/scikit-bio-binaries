/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 2 smoke test: confirms libskbb_wasm.a links and runs under node.
 * Exercises only the non-BLAS surface (get_api_version + permanova_fp64).
 * Full numerical correctness is validated in test_permanova_wasm.cpp
 * (Stage 3) and test_pcoa_wasm.cpp (Stage 6).
 *
 * Exit 0 on success, nonzero on any assertion.
 */

#include "scikit-bio-binaries/util.h"
#include "scikit-bio-binaries/distance.h"

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
    CHECK(skbb_get_api_version() >= 1);

    // Minimal PERMANOVA smoke: 4 samples, 2 groups, no separation.
    // Expected: f-stat finite and non-negative; pvalue in [0, 1].
    const unsigned int n = 4;
    const double mat[] = {
        0.0, 0.2, 0.4, 0.1,
        0.2, 0.0, 0.3, 0.5,
        0.4, 0.3, 0.0, 0.2,
        0.1, 0.5, 0.2, 0.0,
    };
    const unsigned int grouping[] = {0, 0, 1, 1};

    double fstat = -1.0, pvalue = -1.0;
    skbb_permanova_fp64(n, mat, grouping, /*n_perm=*/99, /*seed=*/42,
                        &fstat, &pvalue);

    CHECK(std::isfinite(fstat));
    CHECK(fstat >= 0.0);
    CHECK(pvalue >= 0.0);
    CHECK(pvalue <= 1.0);

    std::printf("smoke: api_version=%u fstat=%.6f pvalue=%.6f failures=%d\n",
                skbb_get_api_version(), fstat, pvalue, failures);
    return failures == 0 ? 0 : 1;
}
