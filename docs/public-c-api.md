# Public C API

The entire shipped surface is three headers, installed to
`$PREFIX/include/scikit-bio-binaries/`, backed by `$PREFIX/lib/libskbb.so`.
Implementations are one-line forwarders in `src/extern/skbb_*.cpp` — never put
logic there.

All matrices are **contiguous row-major `n_dims × n_dims`**, no strides. All
output buffers are **caller-allocated**. Nothing returns an error code; invalid
input is undefined behaviour (see the contracts in
[`permanova.md`](permanova.md#caller-contract-enforced-nowhere--validate-upstream)
and [`pcoa-fsvd.md`](pcoa-fsvd.md#constraints-and-gotchas)).

## Versioning

```c
#define SKBB_API_CURRENT_VERSION 1        /* in util.h — bump on ANY API change */
unsigned int skbb_get_api_version(void);  /* what the loaded .so provides */
```

Each feature group also declares a floor, so a `ctypes`-style consumer can
feature-detect:

| Constant | Covers |
|---|---|
| `SKBB_RANDOM_API_MIN_VERSION` | `skbb_set_random_seed` |
| `SKBB_ACC_API_MIN_VERSION` | the `skbb_*_acc_*` group |
| `SKBB_PERMANOVA_API_MIN_VERSION` | `skbb_permanova_*` |
| `SKBB_CDM_API_MIN_VERSION` | `skbb_center_distance_matrix_*` |
| `SKBB_FSVD_API_MIN_VERSION` | `skbb_fsvd_inplace_*` |
| `SKBB_PCOA_FSVD_API_MIN_VERSION` | `skbb_pcoa_fsvd_*` |

When adding a function: add its `*_MIN_VERSION` constant, bump
`SKBB_API_CURRENT_VERSION`, and add the header to `SHBB_EXTERN_HS` in
`src/Makefile` if the header itself is new. The version check does **not** cover
the dlopen'd GPU plugins.

## `util.h`

```c
unsigned int skbb_get_api_version();
void         skbb_set_random_seed(unsigned int new_seed);   /* global generator */

#define SKBB_ACC_CPU 0
#define SKBB_ACC_NV  1
#define SKBB_ACC_AMD 2
unsigned int skbb_get_acc_mode();                  /* forces detection */
bool         skbb_set_acc_mode(unsigned int t);    /* false if t not compiled in */
bool         skbb_is_acc_reporting();
void         skbb_set_acc_reporting_flag(bool);
```

Under C, `util.h` pulls in `<stdbool.h>`; under C++ it defines `EXTERN extern "C"`.
The x86-level equivalents (`check_use_cpu_x86` etc.) are **internal only** —
`src/util/skbb_detect_acc.hpp`, not exported.

## `distance.h`

```c
void skbb_permanova_fp64(unsigned int n_dims, const double mat[],
                         const unsigned int grouping[], unsigned int n_perm,
                         int seed, double *fstat_ptr, double *pvalue_ptr);
void skbb_permanova_fp32(…float…);
```

* `grouping`: dense 0-based labels, length `n_dims`, `2 <= n_groups < n_dims`.
* `n_perm == 0` → `pvalue = 0.0`.
* `seed < 0` → use the global generator.

## `ordination.h`

```c
/* centering; mat and centered MAY alias */
void skbb_center_distance_matrix_fp64(unsigned int n, const double mat[], double centered[]);
void skbb_center_distance_matrix_fp32(unsigned int n, const float  mat[], float  centered[]);
void skbb_center_distance_matrix_fp64_to_fp32(unsigned int n, const double mat[], float centered[]);

/* eigendecomposition of an ALREADY centered matrix; `centered` is destroyed */
void skbb_fsvd_inplace_fp64(unsigned int n, double centered[], unsigned int n_eighs,
                            int seed, double eigenvalues[], double eigenvectors[]);
void skbb_fsvd_inplace_fp32(…float…);

/* full PCoA */
void skbb_pcoa_fsvd_fp64(unsigned int n, const double mat[], unsigned int n_eighs,
                         int seed, double eigenvalues[], double samples[],
                         double proportion_explained[]);
void skbb_pcoa_fsvd_fp32(…float…);
void skbb_pcoa_fsvd_fp64_to_fp32(unsigned int n, const double mat[], …float outputs…);

/* same, but uses `mat` as the scratch buffer — DESTROYS the input */
void skbb_pcoa_fsvd_inplace_fp64(unsigned int n, double mat[], …);
void skbb_pcoa_fsvd_inplace_fp32(unsigned int n, float  mat[], …);
```

Buffer sizes: `eigenvalues` and `proportion_explained` are `n_eighs`;
`eigenvectors` and `samples` are `n_dims * n_eighs`, **row-major with one row per
sample** (the header's "(n_eighs x n_dims)" wording is misleading; the element
count is the same). `n_eighs <= n_dims`. The three output buffers must not alias
each other — `proportion_explained` is used as scratch during compute.

Eigenvectors are not unit-normalized and their sign is arbitrary.

## Consuming the library

Compiled languages:

```sh
$(CC) $(CFLAGS) my_code.c $(LDFLAGS) -lskbb -o my_exe    # see api_tests/Makefile
```

Python, without any header:

```python
import ctypes
dll = ctypes.CDLL("libskbb.so")
dll.skbb_get_api_version()
```

Static/embedded consumers: `libskbb_inmem.a` (native) or `libskbb_wasm.a`
(emscripten) expose the identical symbol set with no BLAS/LAPACK dependency —
see [`build-system.md`](build-system.md). Downstream emscripten link:

```sh
emcc my_code.c src/libskbb_wasm.a -sEXIT_RUNTIME=1 -I src/extern -o my_code.js
```

## Runtime knobs a caller may care about

`SKBB_USE_GPU`, `SKBB_USE_NVIDIA_GPU`, `SKBB_USE_AMD_GPU`, `SKBB_MAX_CPU`,
`SKBB_GPU_INFO`, `SKBB_CPU_INFO`, `SKBB_TIMING_INFO`, `OMP_NUM_THREADS` —
documented in [`acceleration-dispatch.md`](acceleration-dispatch.md#environment-variables).
`OMP_NUM_THREADS` changes seeded PERMANOVA p-values.
