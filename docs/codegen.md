# Code generation (`src/tools/`)

Every runtime-dispatched kernel exists once in source and is mechanically
expanded into four different shapes by a ~240-line Python line-scanner. This is
the least obvious machinery in the repo; read this before touching
`permanova_dyn_impl.hpp` or `skbb_accapi_impl.hpp`.

## The pieces

```
src/tools/skbb_generate_helper.py     # the engine: print_header() + print_body()
src/tools/generate_permanova_dyn.py   # driver, reads distance/permanova_dyn_impl.hpp
src/tools/generate_skbb_accapi.py     # driver, reads util/skbb_accapi_impl.hpp
```

Invocation (from `src/`, always redirected to a file by the Makefile):

```sh
./tools/generate_permanova_dyn.py <variant> <method> > <out>
#   variant: cpu | cpu_x86_v3 | cpu_x86_v4 | acc_nv | acc_amd   -> namespace skbb_<variant>
#   method:  direct | indirect | api | api_h
```

## The `_T` contract

`print_body()` (`skbb_generate_helper.py:189`) walks the input header line by
line and emits a wrapper for **every function whose declaration line starts with
`static inline ` and contains `_T(`**. Everything else is ignored.

```cpp
static inline int pmn_get_max_parallelism_T() { … }          // -> no-arg wrapper

template<class TFloat>
static inline void pmn_f_stat_sW_T(                          // -> one wrapper per type
        const uint32_t n_dims,
        const TFloat * mat,
        …) { … }
```

Rules the parser enforces (or silently assumes):

* The wrapper name is the text between the third whitespace-separated token and
  `_T(`. The return type is the *second* token (`static inline <type> <name>_T`).
* **A templated (multi-argument) function must return `void`.** Anything else
  raises. Only the zero-argument form may return a value.
* Type expansion is driven by substring match on the argument list:
  `TFloat` → `{float, double}`; `TNum` → `{float, double, uint64_t, uint32_t, bool}`.
  Both are textually replaced by `patch_type()`.
* Argument parsing consumes lines until it sees `") "` — i.e. **the closing paren
  of the parameter list must be followed by a space** (`…, TFloat *out) {`).
  A `)` at end of line, or `){`, breaks the scan.
* The argument *name* is `split()[-1].split('*')[-1]`, so
  `TNum** buf_device` works but a name-less parameter or an array-suffix
  declarator (`TFloat out[]`) will not.

Follow the existing formatting exactly. The parser is positional, not a C++
parser, and it fails by emitting subtly wrong code rather than erroring.

## The four methods

Given `void foo_T(const uint32_t n, TFloat *x)` and variant `acc_nv`
(namespace `skbb_acc_nv`):

| method | Emits | Used for |
|---|---|---|
| `direct` | `template<> void skbb_cpu::foo(uint32_t n, float *x) { foo_T(n, x); }` — plus `#include`s of both the declaration header and the impl header. | CPU variants, where the kernel is compiled by the same compiler into the same library. `permanova_cpu.cpp`, `permanova_cpu_x86_v{3,4}.cpp`, `skbb_accapi_cpu.cpp`. |
| `api_h` | `extern "C" void skbb_acc_nv_foo_float(uint32_t n, float *x);` | The flat C header shared by the plugin and its caller (`permanova_dyn_acc_nv.h`). Emits `#include <stdbool.h>` / `<stdint.h>` first. |
| `api` | `void skbb_acc_nv_foo_float(uint32_t n, float *x) { foo_T(n, x); }` | The **plugin body** — compiled by `nvcc`/`hipcc`/`nvc++`/`amdclang++` into `libskbb_acc_nv.so`. |
| `indirect` | a static `dl_…` function pointer, plus `template<> void skbb_acc_nv::foo(…) { cond_dl_load("skbb_acc_nv_foo_float", &dl_…); (*dl_…)(n, x); }` | The **stub** linked into `libskbb.so`. Also emits `static const char *dl_get_lib_name() { return "libskbb_acc_nv.so"; }` and `#include "util/skbb_dl.cpp"`. |

Special case worth knowing: in `indirect` mode a `bool`-returning no-arg function
whose name contains `found_gpu` gets an extra first line
(`skbb_generate_helper.py:57-59`):

```cpp
if (!dl_load_check()) return false; /* shlib not found */
```

That single hook is what makes a GPU-enabled `libskbb.so` degrade gracefully to
CPU when the plugin `.so` is missing. Every other `indirect` function calls
`cond_dl_load()`, which `exit(1)`s if `dlopen`/`dlsym` fails — so `acc_found_gpu()`
must always be the first plugin call made.

## Adding a dispatched function

1. Write `static inline … myfunc_T(…)` in `permanova_dyn_impl.hpp` (or
   `skbb_accapi_impl.hpp`), guarded by the `#if defined(SKBB_CPU) / SKBB_CUDA /
   SKBB_HIP / OMPGPU / _OPENACC` ladder as appropriate.
2. Declare `myfunc` in the matching `*_dyn.hpp` / `skbb_accapi.hpp` inside
   `namespace SKBB_ACC_NM { … }` — template declaration only, **no definition**.
3. Nothing else. The Makefile rules already regenerate all variants because the
   generated `.cpp`/`.h` depend on the impl header.
4. Call it from the orchestration layer through the right namespace, guarded by
   the runtime dispatch (`src/distance/permanova.cpp:171-225` is the model).

If the new function is called on a device pointer, remember that on CUDA/HIP the
`acc_*_buf` primitives return a *distinct device pointer*, while on
OpenACC/OpenMP-target they return the host pointer with a mapping attached. Code
must work with both — see [`acceleration-dispatch.md`](acceleration-dispatch.md).

## Debugging generated code

The generated files are ordinary text in `src/`. To inspect without a full build:

```sh
cd src && ./tools/generate_permanova_dyn.py acc_nv indirect | less
```

They carry a `// Generated from … // Do not edit by hand` banner. Editing them
is pointless — the next `make` regenerates them, and `make clean` deletes them.
