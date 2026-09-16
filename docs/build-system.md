# Build system

Plain GNU Make, no CMake, no autotools. Two Makefiles plus two included
fragments:

```
Makefile                     # thin dispatcher, `cd src && $(MAKE) <target>`
src/Makefile                 # all native rules; includes the two fragments
src/wasm/emscripten_build.mk # WASM flavor
src/inmem_build.mk           # native-Eigen static-archive flavor
api_tests/Makefile           # black-box tests against the installed .so
api_tests/wasm/Makefile      # same tests under emcc
```

Both fragments are `include`d from `src/Makefile` and **assume `src/` is the
working directory** — every path in them is relative to `src/`.

## Targets

### Top level

| Target | Does |
|---|---|
| `all` | `api` + `install` + `test_bins` |
| `api` | build `libskbb.so` (+ GPU plugins if enabled) |
| `install` | copy shared libs to `$PREFIX/lib`, public headers to `$PREFIX/include/scikit-bio-binaries` |
| `test_bins` | build `src/*.exe` and `api_tests/*.exe` |
| `test` | build + run both test suites |
| `clean`, `clean_install` | see gotchas below |
| `wasm` | `scripts/fetch_eigen.sh` then `src/libskbb_wasm.a` |
| `wasm_test` | build + run the four WASM unit tests under node |
| `wasm_api_test` | build + run `api_tests/` under emcc/node |
| `wasm_clean` | remove WASM artifacts |

`.NOTPARALLEL:` is set at the top level on purpose: the native and WASM DAGs
share generated sources (`permanova_cpu.cpp`, `skbb_accapi_cpu.cpp`), and
`make -j all wasm` would otherwise run the same Python generator twice
concurrently and truncate its own output file. Sub-makes keep their internal
parallelism.

### `src/` only (no top-level passthrough)

| Target | Does |
|---|---|
| `inmem_static` | build `libskbb_inmem.a` |
| `install_inmem` | install the archive + headers |
| `inmem_clean` | remove `libskbb_inmem.a` and `*.inmem.o` |

There is **no `inmem*` target in the root `Makefile`** — invoke it as
`make -C src inmem_static`. (The inmem flavor is newer than the rest of the
build; `src/inmem_build.mk` is a recent addition.)

## Key variables

| Variable | Default | Notes |
|---|---|---|
| `PREFIX` | `$CONDA_PREFIX` | Install root. **If both are unset, `src/Makefile`'s `install` does `mkdir -p /lib`** — no guard. `install_inmem` in `inmem_build.mk` *does* guard. |
| `BLASLIB` | `-llapacke -lcblas` | Conditional (`?=`) so env/CLI can override. Ubuntu needs `-llapacke -lopenblas`; `libcblas.so` under that exact name comes from ATLAS. CI overrides it. |
| `NOGPU` | set automatically on Darwin | Disables `SKBB_ENABLE_ACC_NV/AMD` and the `-ldl` link. |
| `NV_CUDA=Y` | unset | Build the NVIDIA plugin with `nvcc` (`-DSKBB_CUDA`, `-arch=all-major`). |
| `AMD_HIP=Y` | unset | Build the AMD plugin with `hipcc` (`-DSKBB_HIP`, fixed `--offload-arch` list). |
| `NV_CXX` / `AMD_CXX` | unset | Experimental OpenACC (`nvc++`) / OpenMP-target (`amdclang++`) plugins instead. `scripts/enable_{nv,amd}_compiler.sh` set these. |
| `OPT` | empty | Extra flags appended to every compile line. |

Native `CXXFLAGS` are fixed at `-O3 -ffast-math -fopenmp -Wall -std=c++17
-pedantic -I. -fPIC`. Note `-ffast-math`: reassociation is permitted, so
bit-exactness across compilers/optimization levels is not guaranteed for the
float paths — see [`determinism-and-numerics.md`](determinism-and-numerics.md).

On x86_64 non-Darwin, `SKBB_ENABLE_CPU_X86_LEVELS`, `…V3`, `…V4` are defined and
two extra PERMANOVA objects are built with `-march=x86-64-v3` (AVX2) and
`-march=x86-64-v4` (AVX512).

## The `skbb_cpu_tu` macro

`src/Makefile:139-146` defines one canned recipe that emits **three** rules per
translation unit:

```make
$(eval $(call skbb_cpu_tu,<stem>,<source>,<header deps>,<extra flags>))
#   -> <stem>.o        via $(CXX)       $(CXXFLAGS)
#   -> <stem>.wasm.o   via $(WASM_CXX)  $(WASM_CXXFLAGS)
#   -> <stem>.inmem.o  via $(INMEM_CXX) $(INMEM_CXXFLAGS)
```

Use it for any source that is identical across the three flavors. Asymmetric
rules stay explicit: the x86 dispatch objects and the GPU objects are
native-only; `ordination/linalg_backend_lapacke.cpp` (native) vs
`linalg_backend_eigen.cpp` (wasm + inmem) are different sources under different
stems.

After adding a `skbb_cpu_tu` line you must also append the object to the three
object lists — `SKBB_OBJS` (`src/Makefile`), `WASM_OBJS`
(`src/wasm/emscripten_build.mk`), `INMEM_OBJS` (`src/inmem_build.mk`). They are
maintained by hand and nothing checks that they agree.

## Generated sources

Produced by the Python generators at build time; **not committed**; wiped by
`make -C src clean`. See [`codegen.md`](codegen.md).

| Generated file | Recipe |
|---|---|
| `src/permanova_cpu.cpp` | `generate_permanova_dyn.py cpu direct` |
| `src/permanova_cpu_x86_v{3,4}.cpp` | `generate_permanova_dyn.py cpu_x86_v{3,4} direct` |
| `src/permanova_acc_{nv,amd}.cpp` | `… acc_{nv,amd} indirect` (dlopen stubs) |
| `src/permanova_dyn_acc_{nv,amd}.h` | `… acc_{nv,amd} api_h` |
| `src/permanova_dyn_acc_nv.{cu,cpp}` / `…amd.{hip,cpp}` | `… acc_* api` (the plugin body) |
| `src/skbb_accapi_cpu.cpp` and the `skbb_accapi_*` family | `generate_skbb_accapi.py` with the same four methods |

`src/tests/wasm/expected/{permanova,pcoa}_expected.h` are also generated — by
*native* binaries — and are `.gitignore`d on purpose.

## Build outputs

```
src/libskbb_cpu.a      # static, all internal objects; tests link this
src/libskbb.so         # public shared library = extern objects + libskbb_cpu.a + BLAS
src/libskbb_acc_nv.so  # optional NVIDIA plugin, loaded via dlopen at runtime
src/libskbb_acc_amd.so # optional AMD plugin
src/libskbb_wasm.a     # WASM flavor
src/libskbb_inmem.a    # native-Eigen flavor
src/test_{pcoa,permanova}.exe
api_tests/test_{distance,ordination}.exe
```

Note the GPU plugins are **not linked** into `libskbb.so`. `libskbb.so` contains
stub functions that `dlopen("libskbb_acc_nv.so")` on first use; if the file is
absent, `acc_found_gpu()` returns false and the CPU path is taken. That is why
a GPU-capable build still runs on a GPU-less machine.

## WASM specifics

* `scripts/fetch_eigen.sh` downloads Eigen **3.4.0** (SHA-pinned) into
  `.wasm-cache/eigen/include/`. Idempotent; `FORCE_REFRESH=1` re-fetches.
  `scripts/test_eigen_wasm.sh` validates the cache.
* Flags: `-DSKBB_WASM=1 -DSKBB_BLAS_BACKEND_EIGEN=1 -DNOGPU=1 -fno-exceptions
  -Wno-unknown-pragmas`. No `-fopenmp`, no `-pthread`.
  `-Wno-unknown-pragmas` is what lets the `#pragma omp` lines in shared sources
  compile away silently.
* Public headers refer to themselves as `scikit-bio-binaries/<x>.h`; both WASM
  test Makefiles stage copies into a temporary include tree
  (`.wasm-test-include/`, `api_tests/wasm/include/`) to satisfy that.
* Requires `emcc/em++/emar` on `PATH` (emsdk ≥ 5.0.3) and node ≥ 18 to run.

## `clean` gotchas

* `make -C src clean` runs `rm -f *.cpp *.hpp *.h *.cu` **in `src/`**. Safe only
  because every such file there is generated. Do not add hand-written sources at
  that level.
* It does **not** remove `*.hip` (the HIP plugin body), nor the WASM/inmem
  outputs (`wasm_clean` / `inmem_clean` handle those), nor
  `src/tests/wasm/expected/`.
* `clean_install` removes the listed libs and headers from `$PREFIX` — same
  unset-`PREFIX` hazard as `install`.

## CI (`.github/workflows/main.yml`)

Two jobs:

* **build-and-test** — matrix `ubuntu-latest`, `macos-latest`, `ubuntu-24.04-arm`.
  Conda-provisioned compilers + `libcblas liblapacke blas-devel`. On `linux-64`
  only, installs `cuda-compiler` and sets `NV_CUDA=Y` so the CUDA plugin is
  *compiled* (never executed — no GPU on the runner). Tests run with
  `OMP_NUM_THREADS=3` ("a weird number to potentially catch bugs").
* **build-and-test-wasm** — emsdk 5.0.3 + node 20 + a *native* gcc/BLAS toolchain
  (the WASM tests compare against expected values produced by a native binary).
  Runs with `NOGPU=1 BLASLIB="-llapacke -lopenblas"`.

`macos-13` was dropped from the matrix (runners retired).
