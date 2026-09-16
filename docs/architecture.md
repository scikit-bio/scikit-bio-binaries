# Architecture

`scikit-bio-binaries` is a small, self-contained numerical library. It exports a
flat **C ABI** (`skbb_*`) from a shared library `libskbb.so`, implemented in
C++17, with optional per-CPU-generation and per-GPU-vendor variants of the hot
kernels selected **at runtime**.

There are exactly two algorithms today: **PERMANOVA** (`src/distance/`) and
**PCoA via FSVD** (`src/ordination/`). Everything else in the tree exists to
build, dispatch, seed, or test those two.

---

## The layer cake

```
    caller (C, ctypes, emscripten, …)
        │
   ┌────▼──────────────────────────────────────────────┐
   │ 1. Public C ABI          src/extern/*.h + *.cpp   │  extern "C" skbb_*
   │    thin, no logic, just unwraps pointers          │
   └────┬──────────────────────────────────────────────┘
        │
   ┌────▼──────────────────────────────────────────────┐
   │ 2. Algorithm core        skbb::                   │  templated on TFloat,
   │    src/distance/permanova.cpp                     │  instantiated for
   │    src/ordination/principal_coordinate_analysis.cpp│ float + double
   │    owns memory, blocking, RNG, orchestration      │
   └────┬───────────────────────────┬──────────────────┘
        │                           │
   ┌────▼─────────────────────┐ ┌───▼───────────────────┐
   │ 3a. Kernel dispatch      │ │ 3b. Linalg backend    │
   │  skbb_cpu::              │ │  skbb::linalg::       │
   │  skbb_cpu_x86_v3::       │ │  compile-time choice: │
   │  skbb_cpu_x86_v4::       │ │  LAPACKE  |  Eigen    │
   │  skbb_acc_nv::           │ └───────────────────────┘
   │  skbb_acc_amd::          │
   │  chosen at RUNTIME       │
   └────┬─────────────────────┘
        │
   ┌────▼──────────────────────────────────────────────┐
   │ 4. Kernels    distance/permanova_dyn_impl.hpp     │  one source, compiled
   │               util/skbb_accapi_impl.hpp           │  N times with different
   │               #if on SKBB_CPU/CUDA/HIP/OMPGPU/ACC │  flags & compilers
   └───────────────────────────────────────────────────┘
```

Layer 3a is the unusual part and is described in
[The `SKBB_ACC_NM` idiom](#the-skbb_acc_nm-idiom) below and in
[`acceleration-dispatch.md`](acceleration-dispatch.md).

---

## Directory map

| Path | Contents |
|---|---|
| `src/extern/` | Public C headers (`util.h`, `distance.h`, `ordination.h`) + their thin `extern "C"` implementations. **This is the shipped API surface.** Installed to `$PREFIX/include/scikit-bio-binaries/`. |
| `src/distance/` | PERMANOVA. `permanova.{hpp,cpp}` = orchestration; `permanova_dyn.hpp` = the per-variant namespace declaration; `permanova_dyn_impl.hpp` = the actual kernels. |
| `src/ordination/` | PCoA/FSVD. `principal_coordinate_analysis.{hpp,cpp}` + a 3-function linear-algebra backend (`linalg_backend.hpp` with `_lapacke.cpp` / `_eigen.cpp` implementations). |
| `src/util/` | `rand.*` (global + per-thread RNG), `portable_shuffle.hpp`, `skbb_detect_acc.*` (runtime selection), `skbb_accapi*.hpp` (GPU buffer primitives), `skbb_dl.cpp` (dlopen shim), `skbb_dgb_info.hpp` (timing macros). |
| `src/tools/` | The Python code generators. See [`codegen.md`](codegen.md). |
| `src/wasm/`, `src/inmem_build.mk` | Alternate build flavors (Makefile fragments). |
| `src/tests/` | Native C++ tests (hand-rolled harness) + `tests/wasm/` for the WASM parity suite. |
| `api_tests/` | Black-box tests of the installed shared library through the public C headers, compiled as **C99**. |
| `scripts/` | `fetch_eigen.sh` (pinned Eigen drop for WASM/inmem), compiler-enable helpers. |

There are **no sources checked in at the top of `src/`** — everything lives in a
subdirectory. `make -C src clean` runs `rm -f *.cpp *.hpp *.h *.cu` in `src/`
precisely because every file matching those globs there is generated. Never put
a hand-written source directly in `src/`.

---

## The `SKBB_ACC_NM` idiom

The same kernel source is compiled several times, once per hardware variant,
each time into a *different C++ namespace*. The namespace name is injected by a
macro on the compiler command line:

```make
$(CXX) $(CXXFLAGS) -DSKBB_ACC_NM=skbb_cpu_x86_v3 $(X86V3FLAGS) -c permanova_cpu_x86_v3.cpp
```

`distance/permanova_dyn.hpp` and `util/skbb_accapi.hpp` are written to be
included **multiple times in one translation unit**, each time with a different
`SKBB_ACC_NM`:

```cpp
// src/distance/permanova.cpp:19-47
#define SKBB_ACC_NM  skbb_cpu
#include "distance/permanova_dyn.hpp"
#undef SKBB_ACC_NM

#ifdef SKBB_ENABLE_CPU_X86V3
#define SKBB_ACC_NM  skbb_cpu_x86_v3
#include "distance/permanova_dyn.hpp"
#undef SKBB_ACC_NM
#endif
…
```

Two consequences you must not break:

1. **These headers deliberately have no `#ifndef`/`#define` include guard.**
   `permanova_dyn.hpp` and `skbb_accapi.hpp` guard on `#ifdef SKBB_ACC_NM`
   instead (expanding to nothing when the macro is absent). Adding a
   conventional include guard silently reduces them to a single namespace and
   the build fails at link time with missing `skbb_cpu_x86_v3::…` symbols.
2. The declarations are **templates without definitions**; the definitions are
   explicit specializations emitted by the generators into
   `permanova_cpu*.cpp` / `skbb_accapi_*.cpp`. Adding a new type to a kernel
   means teaching the generator about it, not just adding a call site.

## Naming conventions

| Pattern | Meaning |
|---|---|
| `skbb_foo_fp64` / `_fp32` (extern "C") | Public ABI. `fp64_to_fp32` = double input, float output. |
| `skbb::foo` | Internal C++ API, overloaded on `float`/`double`. |
| `foo_T(...)` / `template<class TFloat> foo_T` | An *implementation* function. **The `_T` suffix is load-bearing**: the generators scan for `static inline … _T(` and emit wrappers for exactly those. See [`codegen.md`](codegen.md). |
| `skbb_cpu::`, `skbb_cpu_x86_v3::`, `skbb_cpu_x86_v4::`, `skbb_acc_nv::`, `skbb_acc_amd::` | Per-variant dispatch namespaces (layer 3a). |
| `skbb_acc_nv_pmn_f_stat_sW_double` | The flat `extern "C"` symbol inside a GPU plugin `.so`, resolved by `dlsym`. Shape: `<namespace>_<function>_<type>`. |
| `TFloat` / `TNum` / `TReal` | Template type parameters. `TFloat` → generator emits `{float,double}`; `TNum` → `{float,double,uint64_t,uint32_t,bool}`. |

## Build flavors

Three archives/libraries are produced from largely the same sources. See
[`build-system.md`](build-system.md) for the mechanics.

| Flavor | Output | Linalg | Threads | GPU | x86 dispatch |
|---|---|---|---|---|---|
| native (default) | `libskbb.so` (+ optional `libskbb_acc_nv.so`, `libskbb_acc_amd.so`) | cblas + LAPACKE | OpenMP | dlopen'd plugins | yes, on x86_64 Linux |
| wasm | `src/libskbb_wasm.a` | Eigen (header-only) | none | no | no |
| inmem | `src/libskbb_inmem.a` | Eigen (header-only) | OpenMP | no | no |

The public `skbb_*` symbol set is **identical** across all three. Callers cannot
tell them apart except through `skbb_get_acc_mode()` (always `SKBB_ACC_CPU` for
wasm/inmem) and performance.

## Adding a new algorithm — checklist

1. Core in `src/<area>/<algo>.{hpp,cpp}`, templated `_T` helper + `float`/`double`
   instantiations in namespace `skbb`.
2. Public header `src/extern/<area>.h`: add a `#define SKBB_<NAME>_API_MIN_VERSION`
   constant and the `extern "C"` declarations; bump `SKBB_API_CURRENT_VERSION`
   in `src/extern/util.h`. Add the header to `SHBB_EXTERN_HS` in `src/Makefile`
   if it is new (it is installed and cleaned by name).
3. Thin wrapper in `src/extern/skbb_<area>.cpp`.
4. Register the translation unit **once** with the `skbb_cpu_tu` macro in
   `src/Makefile` (this emits the native `.o`, `.wasm.o` and `.inmem.o` rules)
   and add the object to `SKBB_OBJS`, `WASM_OBJS`, `INMEM_OBJS`.
5. If it needs a runtime-dispatched kernel, follow [`codegen.md`](codegen.md).
6. Tests: native (`src/tests/`), public-API (`api_tests/`), and WASM parity
   (`src/tests/wasm/`) — see [`testing.md`](testing.md).
