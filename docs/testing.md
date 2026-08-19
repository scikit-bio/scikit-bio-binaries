# Testing

Four suites, three harnesses, no framework.

| Suite | Location | Links against | Run with |
|---|---|---|---|
| Native unit | `src/tests/test_{pcoa,permanova}.cpp` | `libskbb_cpu.a` (internal C++ API, `skbb::`) | `make -C src test` |
| Public C API | `api_tests/test_{distance,ordination}.c` | installed `-lskbb`, **C99** | `make -C api_tests test` |
| WASM unit | `src/tests/wasm/test_*_wasm.cpp` | `libskbb_wasm.a` | `make wasm_test` (node) |
| WASM public API | the same `api_tests/*.c` under `emcc` | `libskbb_wasm.a` | `make wasm_api_test` |

`make test` at the top level runs the first two. There is no `ctest`, no Catch2;
the harness is a set of macros in `src/tests/test_helper.hpp` (`SUITE_START`,
`ASSERT`, `SUITE_END`) adapted from BitArray's test harness.

## Native unit tests

`src/tests/test_permanova.cpp` runs each suite **twice** — once before and once
after forcing `skbb::set_use_acc(ACC_CPU)` and (where compiled in)
`set_use_cpu_x86(CPU_X86_BASE)`. That is how the dispatch paths get exercised on
a machine that has AVX512 or a GPU: first pass uses whatever was detected,
second pass uses the baseline. Keep that structure when adding cases.

Expected values come from **scikit-bio** for the statistic (`exp_stat`) and are
loose for the p-value (`|Δ| < 0.05`) because the RNG differs from scikit-bio's.
Tolerance on `fstat` is `1e-5`.

`src/tests/test_pcoa.cpp` compares against hard-coded expected matrices from
scikit-bio (unweighted UniFrac of `test.biom` and `crawford.biom`). Eigenvector
comparisons use `fabs(fabs(got) - fabs(want))` because sign is arbitrary.

Both `main()`s return `tests_failed ? EXIT_FAILURE : EXIT_SUCCESS`.

## Public C API tests

These are the ABI contract tests: compiled as **C99**, including only the
installed public headers, linking `-lskbb`. They require `make install` first
and `$PREFIX/lib` on the loader path. They also assert
`skbb_get_api_version() == SKBB_API_CURRENT_VERSION`, i.e. that the built header
and the built library agree.

> Their `main()` returns `failed ? 1 : 0` where `failed` is reset to `0` at the
> start of every test function — so a failure in an earlier function is masked by
> a later one that passes. `global_failed` is tracked but never used for the
> exit code. See [`known-issues.md`](known-issues.md#api_tests-exit-code-masks-earlier-failures).

## WASM parity suite

The interesting one. It proves the emscripten build behaves like the native
build, and it is the reason several design choices exist (portable shuffle, the
WASM `pmn_get_max_parallelism` value).

```
src/tests/wasm/
  permanova_inputs.hpp   pcoa_inputs.hpp        # fixed inputs, SHARED by generator + test
  generate_permanova_expected.cpp               # NATIVE binary -> prints a C header
  generate_pcoa_expected.cpp
  expected/*.h                                  # GENERATED, .gitignore'd
  test_smoke.cpp  test_permanova_wasm.cpp  test_center_wasm.cpp  test_pcoa_wasm.cpp
```

Flow: `make wasm_test` builds the native `libskbb_cpu.a`, builds and runs the
generators **with `OMP_NUM_THREADS=1`** (enforced twice — by the Makefile rule
and by an explicit check inside the generator, which exits 2 otherwise), writes
`expected/*.h`, then compiles the WASM tests against those headers and runs them
under node.

Expected values are emitted with `%a` (exact hex float) so the round-trip through
the header is lossless, with a `%.17g` decimal in a comment for readability.

**The expected headers are not committed.** That was explicit maintainer
feedback on PR #12 ("Why do we check in the expected files, if they have been
dynamically generated?"). They are in `.gitignore`; regenerate by deleting
`src/tests/wasm/expected/` and re-running `make wasm_test`.

Tolerances (also documented in `README.rst`):

| Check | Tolerance |
|---|---|
| PERMANOVA `fstat`, `pvalue` | **bit-identical** (`==`) |
| `mat_to_centered` | 1e-6 absolute |
| PCoA eigenvalues, proportion_explained | 1e-6 absolute |
| PCoA sample coordinates | 1e-3 absolute, sign-adjusted per axis |

Every WASM test also runs its call **twice with identical inputs** and asserts
bit-equality, catching RNG state leaking between calls.

The sign-adjustment helper in `test_pcoa_wasm.cpp` handles sign flips only. If a
future case has near-degenerate eigenvalues the whole eigenspace can rotate and
a subspace-distance check would be needed instead — the file says so, and the
current 6×6 cases have a well-separated spectrum.

## Adding tests

* New algorithm → add to all four suites. The WASM inputs header must be shared
  by the generator and the test so they cannot drift; both use the
  `static_assert(sizeof(expected)/sizeof(expected[0]) == kCaseCount)` guard.
* Prefer expected values from scikit-bio where a reference exists; state in a
  comment where the number came from.
* If a test produces an incorrect expected value, **do not change the expected
  value** — find out why.
* Run `make test` with a non-power-of-two `OMP_NUM_THREADS` (CI uses 3) — chunk
  arithmetic bugs hide at thread counts of 1.

## Local reproduction of CI

CI's WASM job runs on bare Ubuntu with apt-provided BLAS, not conda, and needs
`BLASLIB="-llapacke -lopenblas"` plus `libblas-dev` for `cblas.h`. A conda env on
the developer machine silently supplies headers the runner lacks — the class of
failure that produced two consecutive CI-fix commits. When touching the link
line, reproduce in a clean container rather than trusting a local conda build.
