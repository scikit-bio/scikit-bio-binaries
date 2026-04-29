/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 6: PCoA (FSVD) correctness under WASM.
 *
 * Replays each PCoACase through the WASM build of skbb_pcoa_fsvd_fp64 and
 * compares against the native (LAPACK) expected values.
 *
 * Tolerance policy:
 *   - eigenvalues: 1e-6 absolute. These are a deterministic function of
 *     the centered distance matrix's dominant eigenspace and should
 *     agree closely across QR/SVD implementations, even though Eigen's
 *     BDCSVD and LAPACK's dgesvd use different bidiagonalization paths.
 *   - proportion_explained: same tolerance (it's eigenvalues / trace).
 *   - samples (coordinates): sign-adjusted per axis, 1e-3 absolute.
 *     Eigenvectors are unique only up to sign (and, in degenerate
 *     eigenspaces, up to rotation). We compare axes independently and
 *     pick the sign that minimizes per-case L2 distance.
 *
 * If a future test case has near-degenerate eigenvalues the sign-only
 * fixup is insufficient and the test would need a subspace-distance
 * helper. The current test cases (6x6 unweighted-UniFrac) have a
 * well-separated spectrum and don't need that.
 */

#include "scikit-bio-binaries/ordination.h"
#include "scikit-bio-binaries/util.h"

#include "tests/wasm/pcoa_inputs.hpp"
#include "tests/wasm/expected/pcoa_expected.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

static int failures = 0;

static bool almost_equal_vec(const char *label, const double *got,
                             const double *want, unsigned int len,
                             double tol) {
    double max_abs_err = 0.0;
    for (unsigned int i = 0; i < len; ++i) {
        double e = std::fabs(got[i] - want[i]);
        if (e > max_abs_err) max_abs_err = e;
    }
    bool ok = max_abs_err <= tol;
    std::printf("    %s  %s  max|delta|=%.3e tol=%.3e len=%u\n",
                label, ok ? "OK" : "FAIL", max_abs_err, tol, len);
    return ok;
}

// Compare two same-shape coordinate matrices (rows x cols, row-major) where
// each column is an eigenvector axis up to sign. For every axis, try both
// sign assignments and accept if either fits within tol.
static bool almost_equal_samples_up_to_sign(const double *got,
                                            const double *want,
                                            unsigned int rows,
                                            unsigned int cols,
                                            double tol) {
    double worst = 0.0;
    unsigned int worst_col = 0;
    int worst_sign = 0;
    for (unsigned int col = 0; col < cols; ++col) {
        double err_pos = 0.0, err_neg = 0.0;
        for (unsigned int row = 0; row < rows; ++row) {
            double g = got[row * cols + col];
            double w = want[row * cols + col];
            double ep = std::fabs(g - w);
            double en = std::fabs(g + w);
            if (ep > err_pos) err_pos = ep;
            if (en > err_neg) err_neg = en;
        }
        double ax_err  = std::fmin(err_pos, err_neg);
        int    ax_sign = (err_pos <= err_neg) ? +1 : -1;
        if (ax_err > worst) {
            worst      = ax_err;
            worst_col  = col;
            worst_sign = ax_sign;
        }
    }
    bool ok = worst <= tol;
    std::printf("    samples  %s  worst-axis=%u sign=%+d max|delta|=%.3e tol=%.3e\n",
                ok ? "OK" : "FAIL", worst_col, worst_sign, worst, tol);
    return ok;
}

int main() {
    using namespace skbb_wasm_test;

    static_assert(sizeof(kPCoAExpected) / sizeof(kPCoAExpected[0])
                  == kPCoACaseCount,
                  "expected-table length must match input-table length");

    const double eigen_tol  = 1e-6;
    const double prop_tol   = 1e-6;
    const double sample_tol = 1e-3;

    for (unsigned int i = 0; i < kPCoACaseCount; ++i) {
        const PCoACase     &c = kPCoACases[i];
        const PCoAExpected &e = kPCoAExpected[i];

        const unsigned int n = c.n;
        const unsigned int k = c.n_eighs;

        double *eigenvalues = new double[k];
        double *samples     = new double[static_cast<uint64_t>(n) * k];
        double *prop        = new double[k];

        // Seed skbb's global RNG for Halko's Gaussian matrix, matching
        // the native generator.
        skbb_set_random_seed(static_cast<unsigned int>(c.seed));
        skbb_pcoa_fsvd_fp64(n, c.mat, k, c.seed,
                            eigenvalues, samples, prop);

        std::printf("  %s\n", c.name);
        bool ok = true;
        ok &= almost_equal_vec("eigenvalues", eigenvalues, e.eigenvalues, k, eigen_tol);
        ok &= almost_equal_vec("prop_expl  ", prop,        e.proportion_explained, k, prop_tol);
        ok &= almost_equal_samples_up_to_sign(samples, e.samples, n, k, sample_tol);

        // Determinism: same seed, same inputs -> identical output on this build.
        double *eigenvalues2 = new double[k];
        double *samples2     = new double[static_cast<uint64_t>(n) * k];
        double *prop2        = new double[k];
        skbb_set_random_seed(static_cast<unsigned int>(c.seed));
        skbb_pcoa_fsvd_fp64(n, c.mat, k, c.seed,
                            eigenvalues2, samples2, prop2);

        bool det_eigen = true, det_samp = true, det_prop = true;
        for (unsigned int j = 0; j < k; ++j) {
            if (eigenvalues[j] != eigenvalues2[j]) det_eigen = false;
            if (prop[j]        != prop2[j])        det_prop  = false;
        }
        for (uint64_t j = 0; j < static_cast<uint64_t>(n) * k; ++j) {
            if (samples[j] != samples2[j]) det_samp = false;
        }
        std::printf("    determinism  %s  (eigen=%s samples=%s prop=%s)\n",
                    (det_eigen && det_samp && det_prop) ? "OK" : "FAIL",
                    det_eigen ? "OK" : "FAIL",
                    det_samp  ? "OK" : "FAIL",
                    det_prop  ? "OK" : "FAIL");
        if (!(det_eigen && det_samp && det_prop)) ok = false;

        if (!ok) ++failures;

        delete[] eigenvalues;  delete[] samples;  delete[] prop;
        delete[] eigenvalues2; delete[] samples2; delete[] prop2;
    }

    std::printf("%s: %u cases, %d failed\n",
                failures == 0 ? "PASS" : "FAIL",
                kPCoACaseCount, failures);
    return failures == 0 ? 0 : 1;
}
