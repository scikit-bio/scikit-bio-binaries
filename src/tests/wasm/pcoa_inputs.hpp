/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Fixed PCoA inputs for native-vs-WASM correctness comparison.
 *
 * For each case we pin a deterministic seed so FSVD's Halko-randomized
 * initial matrix is reproducible. The distance matrices themselves are
 * checked in so the generator and the WASM test read from the same data.
 *
 * The synthetic matrices are built to have a well-separated eigenvalue
 * spectrum (via a dominant rank-k contribution), so per-axis sample
 * coordinates can be compared up to sign with a tight tolerance instead
 * of collapsing into a subspace-only check.
 */

#ifndef SKBB_WASM_PCOA_INPUTS_HPP
#define SKBB_WASM_PCOA_INPUTS_HPP

#include <stdint.h>

namespace skbb_wasm_test {

struct PCoACase {
    const char   *name;
    unsigned int  n;            // matrix dimension (n_samples)
    unsigned int  n_eighs;      // number of eigenvalues requested
    int           seed;
    const double *mat;          // row-major, n*n doubles
};

// ---------------------------------------------------------------------
// Case P1: uuf_6x6 — same 6x6 matrix as the native test_center_mat test.
// ---------------------------------------------------------------------
static const double pcoa_uuf_6x6[] = {
    0.0000000000, 0.2000000000, 0.5714285714, 0.6000000000, 0.5000000000, 0.2000000000,
    0.2000000000, 0.0000000000, 0.4285714286, 0.6666666667, 0.6000000000, 0.3333333333,
    0.5714285714, 0.4285714286, 0.0000000000, 0.7142857143, 0.8571428571, 0.4285714286,
    0.6000000000, 0.6666666667, 0.7142857143, 0.0000000000, 0.3333333333, 0.4000000000,
    0.5000000000, 0.6000000000, 0.8571428571, 0.3333333333, 0.0000000000, 0.6000000000,
    0.2000000000, 0.3333333333, 0.4285714286, 0.4000000000, 0.6000000000, 0.0000000000,
};

static const PCoACase kPCoACases[] = {
    { "uuf_6x6_k3", 6, 3, 12345, pcoa_uuf_6x6 },
    { "uuf_6x6_k2", 6, 2, 12345, pcoa_uuf_6x6 },
};
static const unsigned int kPCoACaseCount =
    sizeof(kPCoACases) / sizeof(kPCoACases[0]);

}  // namespace skbb_wasm_test

#endif
