/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 3: PERMANOVA correctness under WASM.
 *
 * Replays the same fixed inputs (permanova_inputs.hpp) through the WASM
 * build of skbb_permanova_fp64 and compares against the native expected
 * values captured in expected/permanova_expected.h.
 *
 * Tolerance policy:
 *   - fstat: bit-identical (tolerance 0). The unpermuted pseudo-F is a
 *     deterministic function of the distance matrix and grouping with no
 *     RNG dependence; it MUST match across native and WASM.
 *   - pvalue: stochastic variance tolerance. PERMANOVA's p-value is the
 *     fraction of permuted fstats >= observed fstat. The permutation
 *     sequence relies on std::shuffle -> std::uniform_int_distribution,
 *     which is implementation-defined across C++ standard libraries
 *     (libstdc++ vs libc++ differ). std::mt19937 itself is portable, but
 *     how its output is mapped to a permutation is not, so p-values drift
 *     between a gcc-native expected header and an emcc-built WASM run.
 *     We tolerate up to 2 standard errors of the binomial:
 *         sigma = sqrt(p*(1-p) / (n_perm+1))
 *     which is the irreducible Monte-Carlo noise at the chosen n_perm.
 *   - determinism check: we call WASM permanova twice with the same
 *     inputs and assert bit-identical results. This proves the WASM
 *     build itself is deterministic, even if it diverges from a gcc
 *     native build.
 *
 * Exit 0 on success, nonzero on any mismatch.
 */

#include "scikit-bio-binaries/distance.h"

#include "tests/wasm/permanova_inputs.hpp"
#include "tests/wasm/expected/permanova_expected.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

static double pvalue_tolerance(double p, unsigned int n_perm) {
    // Two binomial standard errors, with a floor of 1/(n_perm+1) to
    // guard against zero-width tolerance when p==0 or p==1.
    const double denom = static_cast<double>(n_perm + 1);
    double sigma = std::sqrt(p * (1.0 - p) / denom);
    double floor = 1.0 / denom;
    double tol = 2.0 * sigma;
    return tol > floor ? tol : floor;
}

int main() {
    using namespace skbb_wasm_test;

    static_assert(sizeof(kPermanovaExpected) / sizeof(kPermanovaExpected[0])
                  == kPermanovaCaseCount,
                  "expected-table length must match input-table length");

    int failures = 0;

    for (unsigned int i = 0; i < kPermanovaCaseCount; ++i) {
        const PermanovaCase     &c = kPermanovaCases[i];
        const PermanovaExpected &e = kPermanovaExpected[i];

        double fstat = 0.0, pvalue = 0.0;
        skbb_permanova_fp64(c.n, c.mat, c.grouping,
                            c.n_perm, c.seed, &fstat, &pvalue);

        // Determinism: second run with identical inputs, same seed.
        double fstat2 = 0.0, pvalue2 = 0.0;
        skbb_permanova_fp64(c.n, c.mat, c.grouping,
                            c.n_perm, c.seed, &fstat2, &pvalue2);

        const bool fstat_ok  = (fstat == e.fstat);
        const double ptol    = pvalue_tolerance(e.pvalue, c.n_perm);
        const bool pvalue_ok = (std::fabs(pvalue - e.pvalue) <= ptol);
        const bool det_ok    = (fstat == fstat2) && (pvalue == pvalue2);

        std::printf("  %-12s  fstat %s (got %.17g, want %.17g)  "
                    "pvalue %s (got %.17g, want %.17g, tol=%.4f)  "
                    "det %s\n",
                    c.name,
                    fstat_ok  ? "OK" : "FAIL", fstat,  e.fstat,
                    pvalue_ok ? "OK" : "FAIL", pvalue, e.pvalue, ptol,
                    det_ok    ? "OK" : "FAIL");

        if (!fstat_ok)  ++failures;
        if (!pvalue_ok) ++failures;
        if (!det_ok)    ++failures;
    }

    std::printf("%s: %u cases, %d failures\n",
                failures == 0 ? "PASS" : "FAIL",
                kPermanovaCaseCount, failures);
    return failures == 0 ? 0 : 1;
}
