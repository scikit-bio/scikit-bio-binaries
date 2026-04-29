/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025--, scikit-bio development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Portable Fisher-Yates shuffle built directly on top of std::mt19937.
 *
 * We cannot use std::shuffle: it delegates to std::uniform_int_distribution,
 * whose exact mapping from RNG output to bounded integer is implementation
 * defined (libstdc++, libc++, and MSVC each make a different choice). That
 * means seeded results from std::shuffle are not reproducible across
 * toolchains — which breaks bit-identical PERMANOVA p-values between the
 * native gcc/libstdc++ build and the emscripten/libc++ WASM build.
 *
 * std::mt19937 itself is fully specified by the C++ standard (period,
 * state, output sequence), so building the bounded reduction ourselves
 * removes the only source of cross-toolchain drift in PERMANOVA.
 *
 * The reduction here is the textbook rejection method:
 *    threshold = (2^32 - bound) % bound
 *    loop: r = rng(); if (r < threshold) repeat; else return r % bound
 * which is unbiased and deterministic across platforms.
 */

#ifndef SKBB_PORTABLE_SHUFFLE_HPP
#define SKBB_PORTABLE_SHUFFLE_HPP

#include <stdint.h>

namespace skbb {

// Uniform integer in [0, bound) using rejection on the raw 32-bit output of
// the supplied RNG. Requires bound >= 1.
template <class RNG>
static inline uint32_t portable_uniform_u32(uint32_t bound, RNG &rng) {
    // mt19937 returns values in [0, 2^32). `(uint32_t)(-bound)` is
    // 2^32 - bound in 32-bit wrap-around arithmetic, and modding that by
    // bound gives the smallest rejection threshold that yields a uniform
    // distribution in [0, bound).
    const uint32_t threshold =
        static_cast<uint32_t>(0u - bound) % bound;
    uint32_t r;
    do {
        r = static_cast<uint32_t>(rng());
    } while (r < threshold);
    return r % bound;
}

// In-place Fisher-Yates shuffle of `first[0..n)` driven by `rng`. Matches
// the element order std::shuffle produces on libstdc++ only by coincidence
// — use this explicitly when you want cross-toolchain reproducibility.
template <class T, class RNG>
static inline void portable_shuffle(T *first, uint32_t n, RNG &rng) {
    for (uint32_t i = n; i > 1; --i) {
        const uint32_t j = portable_uniform_u32(i, rng);
        // swap first[i-1], first[j]
        T tmp = first[i - 1];
        first[i - 1] = first[j];
        first[j]     = tmp;
    }
}

}  // namespace skbb

#endif
