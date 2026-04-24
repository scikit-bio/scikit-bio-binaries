/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Shared fixed inputs for the WASM PERMANOVA correctness test.
 * Consumed by both generate_permanova_expected (native) and
 * test_permanova_wasm (WASM), so the two programs cannot drift.
 *
 * Every case pins a seed so std::mt19937 + std::shuffle produce
 * identical permutation sequences across native and WASM, yielding
 * bit-identical fstat and pvalue.
 */

#ifndef SKBB_WASM_PERMANOVA_INPUTS_HPP
#define SKBB_WASM_PERMANOVA_INPUTS_HPP

#include <stdint.h>

namespace skbb_wasm_test {

struct PermanovaCase {
    const char        *name;
    unsigned int       n;
    const double      *mat;
    const unsigned int *grouping;
    unsigned int       n_perm;
    int                seed;
};

// Case 1: 4x4, two obvious groups; matches the smoke test.
static const double case_smoke_mat[] = {
    0.0, 0.2, 0.4, 0.1,
    0.2, 0.0, 0.3, 0.5,
    0.4, 0.3, 0.0, 0.2,
    0.1, 0.5, 0.2, 0.0,
};
static const unsigned int case_smoke_grp[] = {0, 0, 1, 1};

// Case 2: 6x6 unweighted-unifrac-shaped matrix from the existing native
// test suite (src/tests/test_pcoa.cpp:test_center_mat); three groups.
static const double case_uuf_mat[] = {
    0.0000000000, 0.2000000000, 0.5714285714, 0.6000000000, 0.5000000000, 0.2000000000,
    0.2000000000, 0.0000000000, 0.4285714286, 0.6666666667, 0.6000000000, 0.3333333333,
    0.5714285714, 0.4285714286, 0.0000000000, 0.7142857143, 0.8571428571, 0.4285714286,
    0.6000000000, 0.6666666667, 0.7142857143, 0.0000000000, 0.3333333333, 0.4000000000,
    0.5000000000, 0.6000000000, 0.8571428571, 0.3333333333, 0.0000000000, 0.6000000000,
    0.2000000000, 0.3333333333, 0.4285714286, 0.4000000000, 0.6000000000, 0.0000000000,
};
static const unsigned int case_uuf_grp[] = {0, 0, 1, 1, 2, 2};

// Case 3: 8x8 equal-group no-signal (ties), from tests/test_permanova.cpp style.
static const double case_ties_mat[] = {
    0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5,
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5,
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0,
};
static const unsigned int case_ties_grp[] = {0, 0, 0, 0, 1, 1, 1, 1};

static const PermanovaCase kPermanovaCases[] = {
    { "smoke_4x4",  4, case_smoke_mat, case_smoke_grp,  99, 42 },
    { "uuf_6x6",    6, case_uuf_mat,   case_uuf_grp,   199,  1 },
    { "ties_8x8",   8, case_ties_mat,  case_ties_grp,  499,  7 },
};
static const unsigned int kPermanovaCaseCount =
    sizeof(kPermanovaCases) / sizeof(kPermanovaCases[0]);

}  // namespace skbb_wasm_test

#endif
