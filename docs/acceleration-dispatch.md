# Acceleration & runtime dispatch

Two orthogonal selections happen at runtime, both implemented in
`src/util/skbb_detect_acc.cpp`:

1. **Accelerator**: CPU / NVIDIA / AMD (`check_use_acc()`, line 65).
2. **x86 microarchitecture level**: base / v3 (AVX2) / v4 (AVX512)
   (`check_use_cpu_x86()`, line 214). Compiled in only when
   `SKBB_ENABLE_CPU_X86_LEVELS` is defined — x86_64 Linux only.

Both cache their answer in a file-static `int` initialized to `-1`, so detection
runs once per process. The comment `// we can assume int is atomic` is the
extent of the thread-safety argument: a benign race where two threads both run
detection and store the same value.

## Accelerator selection

```
detected = ACC_CPU
if SKBB_ENABLE_ACC_NV  and skbb_acc_nv::acc_found_gpu():   detected = ACC_NV
    unless SKBB_USE_NVIDIA_GPU ∈ {N,NO,n,no,NEVER,never}
if SKBB_ENABLE_ACC_AMD and skbb_acc_amd::acc_found_gpu():  detected = ACC_AMD   (only if still CPU)
    unless SKBB_USE_AMD_GPU ∈ {…}
if SKBB_USE_GPU ∈ {N,NO,n,no,NEVER,never}:                 detected = ACC_CPU
```

NVIDIA wins over AMD when both are present. `acc_found_gpu()` for a variant is a
`dlopen` of that variant's plugin — see below.

Public overrides (`src/extern/util.h`):

* `skbb_get_acc_mode()` → `SKBB_ACC_CPU|NV|AMD` (0/1/2), forcing detection.
* `skbb_set_acc_mode(t)` → returns `false` and changes nothing if `t` names a
  variant this build wasn't compiled with. Sets the cached value, so it also
  works as a "skip detection" switch.
* `skbb_is_acc_reporting()` / `skbb_set_acc_reporting_flag(bool)`.

## Environment variables

| Variable | Effect |
|---|---|
| `SKBB_USE_GPU=N` | force CPU (any accelerator) |
| `SKBB_USE_NVIDIA_GPU=N` / `SKBB_USE_AMD_GPU=N` | disable one vendor |
| `SKBB_GPU_INFO=Y` | print accelerator selection decisions to stdout; also makes `skbb_dl.cpp` announce which plugin `.so` it loaded |
| `SKBB_MAX_CPU=basic\|base\|x86-64-v2\|sse\|sse3\|sse4\|avx` | cap at base; `x86-64-v3`/`avx2` caps at v3. Unknown strings are ignored silently. |
| `SKBB_CPU_INFO=Y` | print x86 level selection |
| `SKBB_TIMING_INFO=Y` | per-phase wall-clock via the `SETUP_TDBG`/`TDBG_STEP` macros in `util/skbb_dgb_info.hpp` |
| `OMP_NUM_THREADS` | standard OpenMP. **Changes seeded PERMANOVA p-values** — see [`determinism-and-numerics.md`](determinism-and-numerics.md). |

The truthiness convention is inverted and worth noting: for the `*_INFO`
variables, *any* value enables reporting **except** the explicit negatives
(`NO`, `N`, `no`, `n`, `NEVER`, `never`). `SKBB_GPU_INFO=0` therefore turns
reporting **on**.

One inconsistency: `skbb_dl.cpp:31` and `:74` check `env_cpu_info[0]=='Y'`
rather than using the same negative-list convention, so `SKBB_GPU_INFO=y` prints
the selection messages but not the "Using shared library …" line.

## The dlopen plugin mechanism

`libskbb.so` never links CUDA/HIP. Instead, `generate_*.py … indirect` emits
stubs that lazily resolve symbols out of a sibling `.so`:

```
libskbb.so
  skbb_acc_nv::acc_found_gpu()      -> dl_load_check()      -> dlopen("libskbb_acc_nv.so")
  skbb_acc_nv::pmn_f_stat_sW<double>-> cond_dl_load("skbb_acc_nv_pmn_f_stat_sW_double")
                                                            -> dlsym
```

`src/util/skbb_dl.cpp` (textually `#include`d into each generated stub file, so
each has its own `dl_handle` and mutex):

* `dl_load_check()` — `dlopen` with `RTLD_LAZY`; returns `false` if the library
  is absent. **Only `acc_found_gpu` uses this soft path.**
* `cond_dl_load(name, &ptr)` — mutex-guarded lazy `dlsym`; on any failure it
  prints to stderr and `exit(1)`. Hard failure by design: if detection said the
  GPU is there, a missing symbol is a build inconsistency.
* The `.so` is found through the normal loader search path (`LD_LIBRARY_PATH`,
  rpath); there is no absolute path or version suffix — `dl_get_lib_name()`
  returns the bare `"libskbb_acc_nv.so"`.

Consequence: a plugin built against a different `skbb` version will be loaded
happily as long as the symbol names match. The `skbb_get_api_version()`
mechanism does **not** cover the plugin boundary.

## The accelerator buffer API (`skbb_accapi`)

`src/util/skbb_accapi_impl.hpp` is a 6-primitive abstraction over "move a buffer
to the device", specialized by `#if` on `SKBB_CUDA` / `SKBB_HIP` / `OMPGPU` /
`_OPENACC` / nothing:

| Primitive | CUDA / HIP | OpenACC / OpenMP-target | plain CPU |
|---|---|---|---|
| `acc_create_buf(host, &dev, n)` | `cudaMalloc` → **new pointer** | `enter data create`, `dev = host` | `dev = host` |
| `acc_copyin_buf(host, &dev, n)` | `cudaMalloc` + `Memcpy H2D` | `enter data copyin`, `dev = host` | `dev = host` |
| `acc_update_device(host, dev, start, end)` | `Memcpy` of the slice | `update device(host[start:end])` | no-op |
| `acc_copyout_buf(host, dev, n)` | `Memcpy D2H` + free | `exit data copyout` | no-op |
| `acc_destroy_buf(dev, n)` | `cudaFree` | `exit data delete` | no-op |
| `acc_found_gpu()`, `acc_need_alt()`, `acc_wait()` | device count / sync | ditto | `false` / `false` / no-op |

**The critical asymmetry**: under CUDA/HIP the device pointer is a *distinct*
allocation and host memory is untouched until `copyout`; under
OpenACC/OpenMP-target the device pointer *is* the host pointer with a mapping
attached. Orchestration code must therefore keep both `x` and `x_device`
variables and never assume they differ (see
`src/distance/permanova.cpp:123-148`).

CUDA/HIP failures `throw std::runtime_error`. Nothing in `libskbb.so` catches
them, and the public C ABI has no error channel — a CUDA OOM propagates out of
an `extern "C"` function, which is UB for C callers. Treat GPU paths as
abort-on-failure.

## x86 level dispatch

Only PERMANOVA's inner kernel is dispatched this way — three copies of
`pmn_f_stat_sW` are compiled (`-march=x86-64-v3`, `-march=x86-64-v4`, and
baseline) into `skbb_cpu_x86_v3::`, `skbb_cpu_x86_v4::`, `skbb_cpu::`. Selection
uses `__builtin_cpu_init()` + `__builtin_cpu_supports("x86-64-v4"/"v3")`.

`pmn_get_max_parallelism()` is deliberately taken from `skbb_cpu::` regardless of
the selected level (`src/distance/permanova.cpp:97`) — the chunk size must not
depend on the microarchitecture, or seeded results would differ between machines.

PCoA is *not* x86-dispatched; it relies on `#ifdef __AVX2__`/`__AVX__` inside
`E_matrix_means` (`principal_coordinate_analysis.cpp:58-114`), evaluated against
the flags of the single baseline compile — i.e. **the unrolled paths are dead
code in the default build** and BLAS does the heavy lifting anyway.

## What wasm/inmem see

Neither defines `SKBB_ENABLE_ACC_NV/AMD` nor `SKBB_ENABLE_CPU_X86_LEVELS`, so all
the `#if`-guarded blocks vanish: `check_use_acc()` always returns `ACC_CPU`,
`skbb_set_acc_mode(SKBB_ACC_NV)` returns `false`, and the only surviving
`skbb_accapi` implementation is the identity (`dev = host`). The `skbb_accapi_cpu`
objects are still built and linked — they're the no-op path.
