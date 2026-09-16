# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ABSOLUTE HARD REQUIREMENTS
- NEVER use `rm` without permission
- NEVER commit generated files (`src/*.cpp` at that level, `src/tests/wasm/expected/`, `*.o`, `*.a`)
- NEVER `git add -A`, `git add .`, or `git commit -a` — stage explicit paths

## Project Overview

**scikit-bio-binaries** provides optimized C/C++ implementations of numerical
algorithms in support of `scikit-bio`, exposed as a flat C ABI (`skbb_*`) from
`libskbb.so`. Two algorithms today: **PERMANOVA** and **PCoA via FSVD**.

The distinguishing property is that hot kernels are compiled multiple times —
per x86 microarchitecture level and per GPU vendor — and selected **at runtime**,
with GPU support arriving through `dlopen`'d plugin libraries so a GPU-enabled
build still runs on a GPU-less machine. Three build flavors (native shared
library, WebAssembly static archive, native Eigen-only static archive) share one
set of sources and expose an identical symbol set.

This repo is upstream `scikit-bio/scikit-bio-binaries`, maintained by
@sfiligoi. We contribute via PRs from the `the-miint` fork.

## Priorities

1. Red/green/refactor Test Driven Development (TDD)
2. Verifiably correct code
3. Maintainable code, using Don't Repeat Yourself (DRY) and Keep It Simple Stupid (KISS)
4. Performance

## Rules

These rules apply to every task in this project unless explicitly overridden.

### Rule 1 — Think Before Coding
Bias: caution over speed on non-trivial work.
State assumptions explicitly. Ask rather than guess.
Push back when a simpler approach exists. Stop when confused.

### Rule 2 — Simplicity First
Minimum code that solves the problem. Nothing speculative.
No abstractions for single-use code.

### Rule 3 — Surgical Changes
Touch only what you must. Don't improve adjacent code.
Match existing style. Don't refactor what isn't broken.

### Rule 4 — Goal-Driven Execution
Define success criteria. Loop until verified.

### Rule 5 — Surface conflicts, don't average them
If two patterns contradict, pick one (more recent / more tested).
Explain why. Flag the other for cleanup.

### Rule 6 — Read before you write
Before adding code, read exports, immediate callers, shared utilities.
If unsure why existing code is structured a certain way, ask.

### Rule 7 — Tests verify intent, not just behavior
Tests must encode WHY behavior matters, not just WHAT it does.
A test that can't fail when business logic changes is wrong.

### Rule 8 — Checkpoint after every significant step
Summarize what was done, what's verified, what's left.
Don't continue from a state you can't describe back.

### Rule 9 — Match the codebase's conventions, even if you disagree
Conformance > taste inside the codebase.
If you think a convention is harmful, surface it. Don't fork silently.

### Rule 10 — Fail loud
"Completed" is wrong if anything was skipped silently.
Don't silently skip tests you caused to be skipped.
Default to surfacing uncertainty, not hiding it.

## Style, churn, and comments — non-negotiable here

There is **no formatter and no format check in this repo**. Nothing will undo a
gratuitous reformat, and every extra line is reviewer cost paid by the upstream
maintainer. Full detail in **[`docs/code-style.md`](docs/code-style.md)** — read
it before your first edit in a session that writes code. The three rules that
matter most:

1. **Match the file you are editing.** Indentation is genuinely inconsistent
   across this tree (1-space, 2-space, 4-space, tabs) and there is no globally
   correct answer. Copy the local convention.
2. **No unnecessary churn.** No reformatting, re-wrapping, re-ordering,
   typo-fixing, or modernizing outside the functional change you were asked for.
   Unrelated fixes — including anything in
   [`docs/known-issues.md`](docs/known-issues.md) — get their own commit and
   their own test.
3. **Comments must be succinct and durable.** Explain *why*, in the shortest form
   that survives. Write for someone reading the file in two years who knows
   nothing about the current change. **Never narrate the work**: no "Stage 3",
   no "per reviewer guidance", no "previous fix tried…", no PR numbers. That
   belongs in the commit message. The tree already contains such comments and
   they have gone stale and now mislead.

Standing maintainer feedback (from PR review history): **do not duplicate build
infrastructure** — add a translation unit in one place, via the `skbb_cpu_tu`
macro; **do not commit generated artifacts**; and prefer the boring solution,
because complexity in the build/dispatch layer gets flagged.

## Build Commands

Conda-provided compilers recommended; system ones work.

```bash
conda create -n skbb-build -c conda-forge gxx_linux-64      # or clangxx_osx-* on macOS
conda activate skbb-build
conda install -c conda-forge libcblas liblapacke blas-devel make

make clean && make clean_install && make all   # api + install + test binaries
make test                                      # native unit tests + public C API tests
```

WebAssembly (needs an activated emsdk ≥ 5.0.3 and node ≥ 18):

```bash
scripts/fetch_eigen.sh    # pinned Eigen 3.4.0 into .wasm-cache/
make wasm                 # src/libskbb_wasm.a
make wasm_test            # smoke, PERMANOVA, centering, PCoA under node
make wasm_api_test        # public C API parity
```

Native Eigen-only static archive (no top-level target):

```bash
make -C src inmem_static  # src/libskbb_inmem.a
```

GPU builds are off by default: `export NV_CUDA=Y` (needs `cuda-compiler`) or
`export AMD_HIP=Y` (needs `hipcc`).

Useful at runtime: `SKBB_USE_GPU=N`, `SKBB_MAX_CPU=basic`, `SKBB_GPU_INFO=Y`,
`SKBB_CPU_INFO=Y`, `SKBB_TIMING_INFO=Y`, `OMP_NUM_THREADS`.

## Testing

If a test produces an **incorrect expected value**: DO NOT change the expected
value without permission.

Run `make test` with a non-power-of-two `OMP_NUM_THREADS` (CI uses 3) — chunking
bugs hide at 1 thread. `OMP_NUM_THREADS` changes seeded PERMANOVA p-values by
design; see [`docs/determinism-and-numerics.md`](docs/determinism-and-numerics.md).

## Architecture and Patterns (detailed docs)

Load only what the task needs.

- **[`docs/architecture.md`](docs/architecture.md)** — the layer cake (public C ABI → `skbb::` core → per-variant dispatch namespaces → kernels), directory map, naming conventions, the three build flavors, and the **`SKBB_ACC_NM` multi-include idiom** (headers deliberately without include guards — adding one breaks the build in a confusing way). Start here.
- **[`docs/code-style.md`](docs/code-style.md)** — per-file conventions, comment policy, what the maintainer has pushed back on. Read before writing code.
- **[`docs/build-system.md`](docs/build-system.md)** — Make targets and variables, the `skbb_cpu_tu` canned recipe, generated sources, install/clean footguns, CI layout. Read when touching any Makefile or adding a translation unit.
- **[`docs/codegen.md`](docs/codegen.md)** — the Python generators that expand `_T`-suffixed kernels into `direct`/`indirect`/`api`/`api_h` wrappers. Read before editing `permanova_dyn_impl.hpp` or `skbb_accapi_impl.hpp`; the parser is positional and fails silently.
- **[`docs/acceleration-dispatch.md`](docs/acceleration-dispatch.md)** — runtime CPU/GPU selection, all environment variables, the `dlopen` plugin mechanism, and the device-buffer API (CUDA/HIP return a *distinct* pointer; OpenACC/OpenMP-target do not).
- **[`docs/permanova.md`](docs/permanova.md)** — the statistic, the chunked permutation loop, `NBLOCK`/`TILE` blocking, the CUDA and OpenACC kernels, and the caller contract (dense 0-based groupings) that nothing validates.
- **[`docs/pcoa-fsvd.md`](docs/pcoa-fsvd.md)** — centering math, the six FSVD steps, the three-function linalg backend and its **column-major `svd_no` stride contract**, output layouts, and the aliasing rules between output buffers.
- **[`docs/determinism-and-numerics.md`](docs/determinism-and-numerics.md)** — what is bit-reproducible and what is not, why `std::shuffle` is banned, and the list of changes that silently alter seeded p-values. Read before touching RNG, chunk sizes, or accumulation order.
- **[`docs/testing.md`](docs/testing.md)** — the four suites, the native→WASM expected-value generation flow, tolerances, and how to add tests.
- **[`docs/public-c-api.md`](docs/public-c-api.md)** — every exported symbol, buffer sizes and ownership, and the API-versioning rules to follow when adding a function.
- **[`docs/known-issues.md`](docs/known-issues.md)** — verified defects (including a real PERMANOVA correctness bug in the partial-chunk path), suspected ones, and latent build fragilities. Check here before debugging something odd; do not fix these as drive-bys.
