# Determinism and numerics

This library's reproducibility guarantees are narrower than they look. Read this
before changing anything that touches RNG, thread counts, chunk sizes, or the
order of floating-point accumulation.

## The RNG model

`src/util/rand.{hpp,cpp}`. `skbb::TRAND` is `std::mt19937` — chosen because the
C++ standard fully specifies its state and output sequence, so it is identical
across libstdc++, libc++, and emscripten.

Two levels:

* **Global generator** — one file-static `std::mt19937`, seeded by
  `skbb_set_random_seed()`. Used whenever a `seed` argument is `< 0`.
  Process-wide mutable state: not thread-safe, and results depend on how many
  random numbers earlier calls consumed.
* **`RandomGeneratorArray(size, seed)`** — an array of independent generators,
  one per parallel slot. `rga_init` (rand.cpp:28) draws `size` seeds *in order*
  from either a locally-seeded generator (`seed >= 0`) or the global one
  (`seed < 0`), then seeds each element. This is what makes threaded PERMANOVA
  deterministic: slot *i* always gets the same stream regardless of scheduling.

## `portable_shuffle`, and why `std::shuffle` is banned

`src/util/portable_shuffle.hpp`. `std::shuffle` delegates to
`std::uniform_int_distribution`, whose mapping from RNG output to a bounded
integer is **implementation-defined**; libstdc++, libc++ and MSVC each choose
differently. Seeded results would then differ between the native gcc build and
the emscripten build.

The replacement is a textbook rejection-method Fisher-Yates:

```
threshold = (2^32 - bound) % bound
do { r = rng(); } while (r < threshold);
return r % bound;
```

Unbiased and byte-identical everywhere. **Do not reintroduce `std::shuffle`,
`std::uniform_int_distribution`, or any `<random>` distribution on a hot,
reproducibility-sensitive path.** (`std::normal_distribution` *is* used in
FSVD's Halko projection — see the caveat below.)

## What is guaranteed, and what is not

| Property | Guaranteed? |
|---|---|
| Same build, same inputs, same seed, same thread count → identical output | **Yes.** The WASM tests assert bit-equality on a repeat call. |
| Native vs WASM PERMANOVA `fstat`/`pvalue` at a fixed seed | **Yes, bit-identical** — same arithmetic, same mt19937, portable shuffle, and `pmn_get_max_parallelism()` returns 32 under WASM to match native single-thread. |
| Same seed across **different `OMP_NUM_THREADS`** | **No.** `PERM_CHUNK = 2·threads·16` changes `step_perms`, which changes how many generators are created and how permutations are chunked. |
| Native vs WASM PCoA | Approximately: 1e-6 on eigenvalues/proportions, 1e-3 on coordinates (sign-adjusted). LAPACK `dgesvd` and Eigen `BDCSVD` bidiagonalize differently. Observed drift in the shipped cases is ~1e-15; the headroom is for larger/ill-conditioned inputs. |
| FSVD seeded results across toolchains | **No.** `std::normal_distribution` has the same implementation-defined problem as `uniform_int_distribution`. The WASM PCoA test therefore compares with tolerances, not bit-equality. |
| fp32 vs fp64 agreement | Only to fp32 precision. Note accumulators are `double` in both PERMANOVA kernels regardless of `TFloat`. |
| Bit-stability across compiler/optimization changes | **No.** `-ffast-math` is on for the native build, permitting reassociation. |

## Things that silently change seeded results

Treat any change to these as an API-visible change and regenerate the WASM
expected values (see [`testing.md`](testing.md)):

1. `pmn_get_max_parallelism()` for any backend.
2. `NBLOCK` (16) or the tail decomposition in `pmn_f_stat_sW_cpu`.
3. The number, order, or seeding of `RandomGeneratorArray` elements.
4. Where `portable_shuffle` is called, or how many draws it makes.
5. Choosing a different x86 level for `pmn_get_max_parallelism` (currently
   deliberately always `skbb_cpu::`, permanova.cpp:97).
6. The order of the permutation loop / whether rows are reset between chunks.

The `#error` in `pmn_get_max_parallelism_T` for non-OpenMP non-WASM CPU builds
exists for exactly this reason: an arbitrary fallback chunk size would silently
change every seeded p-value.

## Precision choices in the kernels

* PERMANOVA accumulates `s_W` in `double` even for the `float` instantiation
  ("Use full precision for intermediate compute, to minimize accumulation
  errors"), at both the per-row and per-tile level, and in the CUDA shared-memory
  reduction. `TFloat` controls only the storage type of the distance matrix.
* `sum_upper_square` likewise reduces into a `double` and casts at the end.
* PCoA centering accumulates row sums in `TReal` (the *output* type), so the
  fp64→fp32 mixed path accumulates in fp32.
* The CPU and GPU PERMANOVA kernels apply `inv_group_sizes` at different
  granularity (per row vs per element), so their results differ in the last bits.

## Seed semantics at the public API

Every seeded entry point takes `int seed`:

* `seed >= 0` → deterministic, self-contained: a local generator is seeded with
  it and nothing touches global state.
* `seed < 0` → draw from the library-global generator (settable with
  `skbb_set_random_seed`). Reproducible only if you control every prior call in
  the process.

Note `skbb_set_random_seed` takes `unsigned int` while the per-call `seed` is
`int`, so the per-call seed space is half as large.
