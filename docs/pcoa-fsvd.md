# PCoA and FSVD

Files: `src/ordination/principal_coordinate_analysis.{hpp,cpp}`,
`linalg_backend.hpp` + `linalg_backend_{lapacke,eigen}.cpp`.
Public entries: `skbb_center_distance_matrix_*`, `skbb_fsvd_inplace_*`,
`skbb_pcoa_fsvd_*`.

## Pipeline

```
mat (n×n distances)
  │  mat_to_centered            Legendre & Legendre 1998, eq. 9.20 + 9.21
  ▼
centered (n×n)  ── trace kept for proportion_explained
  │  find_eigens_fast           Halko et al. 2011 randomized SVD ("FSVD")
  ▼
eigenvalues[k], eigenvectors (n×k)
  │  scale each axis by sqrt(eigenvalue)
  ▼
samples (n×k), proportion_explained[k] = eigenvalues / trace(centered)
```

## Centering (`mat_to_centered`)

Two passes, both OpenMP-parallel over rows:

1. `E_matrix_means` (line 34) — `E[i,j] = -0.5 · D[i,j]²`, accumulating row sums
   and the global sum in the same pass. `row_means[i] = rowsum/n`,
   `global_mean = (globalsum/n)/n`.
2. `F_matrix_inplace` (line 133) — `F[i,j] = E[i,j] - row_mean[i] - row_mean[j]
   + global_mean`, tiled 512×512 to keep `row_means` hot.

Since the matrix is symmetric, column means are row means — that identity is why
only one mean vector is computed.

`mat` and `centered` may alias (in-place). The `E_matrix_means` `#ifdef __AVX2__`
/ `__AVX__` 8×/4× unrolled variants are compile-time only and, with the default
flags (no `-march`), the scalar tail loop is what actually runs.

The mixed overload `mat_to_centered(n, const double*, float*)` reads fp64 and
writes fp32, computing means in fp32 (`TReal` = the *output* type). That is the
`skbb_center_distance_matrix_fp64_to_fp32` path.

## FSVD (`find_eigens_fast_T`, line 360)

Randomized eigendecomposition of the (symmetric, PSD-ish) centered matrix.
With `k = n_eighs + 2`:

| Step | Code | Shapes |
|---|---|---|
| 1. Random projection with one power iteration | `centered_randomize_T` (line 211) | `G` ~ N(0,1), n×k. `H = [A·G ; A·A·(A·G)]`, n×2k |
| 2. Orthonormalize | `QR` class (line 254) → `qr_inplace` | `Q` = n×qcols, `qcols = min(n, 2k)`; H buffer becomes Q |
| 3. Project | `qdot_r_sq` → `gemm_nn` | `T = A·Q`, n×qcols (uses `Aᵀ == A`) |
| 4. SVD | `svd_it_T` → `svd_no` | singular values → `S`; `Vᵀ` overwrites `T` |
| 5. Lift back | `qdot_l_sq` | `U = Q·V`, n×qcols, reusing the `T` buffer |
| 6. Truncate + transpose | `transpose_T` | first `n_eighs` values/columns → outputs |

The power iteration is fixed at **one** pass (`A·A·(A·G)` — three gemms). There
is no convergence check and no oversampling knob beyond the hard-coded `+2`.

Two ownership subtleties: `QR` takes ownership of the `H` malloc and `free`s it
in its destructor; `Ut` then takes ownership of the `T` malloc. Both are raw
`malloc`/`free` interleaved with `new[]`/`delete[]` elsewhere in the file — match
whatever the surrounding function already uses.

`qr_i_T` and `svd_it_T` failures call `exit(1)` after printing to stderr
("should never fail"). There is no error return to the caller.

## Linear-algebra backend

`skbb::linalg` is deliberately three functions, not a BLAS facade
(`linalg_backend.hpp`):

```cpp
void gemm_nn(m, n, k, A, B, C);         // C = A·B, column-major, alpha=1 beta=0
int  qr_inplace(rows, cols, H, &qcols); // H -> thin Q (rows × min(rows,cols))
int  svd_no(rows, cols, T, S);          // LAPACK jobu='N', jobvt='O'
```

Compile-time selection: `-DSKBB_BLAS_BACKEND_EIGEN=1` picks
`linalg_backend_eigen.cpp` (WASM + inmem); otherwise the native build compiles
`linalg_backend_lapacke.cpp` (`cblas_?gemm`, `LAPACKE_?geqrf`+`?orgqr`,
`LAPACKE_?gesvd`). Both are compiled under the object stem
`ord_linalg_backend*.o`; only one is linked into a given artifact.

**The `svd_no` output layout contract is the sharp edge.** `Vᵀ` occupies the
first `k = min(rows,cols)` rows of the `rows × cols` column-major buffer, so it
must be read with stride `rows`, not packed. Rows `[k, rows)` are undefined.
`transpose_sq_st_T(qr_obj.cols, qr_obj.rows, T, W)` exists exactly to honour
that stride. The Eigen implementation reproduces the layout by hand
(`linalg_backend_eigen.cpp:104-110`) rather than inheriting it.

Everything here is **column-major** ("FORTRAN-style ColOrder"), while the
public matrices and the final outputs are row-major. The transposes at the
boundary are not redundant.

## Output layouts and sizes

| Buffer | Size | Layout |
|---|---|---|
| `mat` / `centered` | `n*n` | row-major (symmetric, so also column-major) |
| `eigenvalues` | `n_eighs` | descending |
| `eigenvectors` (from `skbb_fsvd_inplace_*`) | `n * n_eighs` | **row-major, n rows × n_eighs cols** |
| `samples` (from `skbb_pcoa_fsvd_*`) | `n * n_eighs` | **row-major, one row per sample** |
| `proportion_explained` | `n_eighs` | `eigenvalue / trace(centered)` |

The header comment says "Matrix of size (n_eighs x n_dims)" while the prose says
"row-indexed by the sample id"; the code (`transpose_T(n_dims, n_eighs, …)` and
the scaling loop at line 575) is unambiguous: **n rows of n_eighs values**. Total
element count is the same either way, which is why the tests never caught the
wording.

## Constraints and gotchas

* `n_eighs <= n_dims` is required (`eigenvalues[i] = S[i]` for `i < n_eighs`,
  and `S` only holds `min(n, 2k)` meaningful values). Not checked.
* **`proportion_explained` is used as scratch** inside `pcoa_T` (line 567) to
  hold `sqrt(eigenvalues)` before being overwritten with the real values. The
  three output buffers must not alias each other.
* `samples` **is** the eigenvector buffer during compute (line 531) — that is
  why it needs `n*n_eighs` elements even before scaling.
* `skbb_fsvd_inplace_*` and `skbb_pcoa_fsvd_inplace_*` destroy their input
  buffer. The `NewCentered` / `InPlaceCentered` pair (lines 441-484) is the only
  difference between the two families: one allocates an `n*n` scratch, the other
  reuses `mat`.
* Eigenvectors are **not normalized to unit length** on output (documented in
  the public header) and their sign is arbitrary — compare up to sign.
* Seeding: `seed < 0` means "draw from the library-global generator", which is
  process-shared mutable state. See
  [`determinism-and-numerics.md`](determinism-and-numerics.md).
