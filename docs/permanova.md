# PERMANOVA

Files: `src/distance/permanova.{hpp,cpp}`, `permanova_dyn.hpp`,
`permanova_dyn_impl.hpp`; public entry `skbb_permanova_fp{64,32}`.

## The statistic

Given a symmetric distance matrix `D` (n×n) and a grouping vector `g` (length n,
values `0..n_groups-1`):

```
s_T   = ( Σ_{i<j} D[i,j]² ) / n                       # total sum of squares
s_W   = Σ_{i<j, g[i]==g[j]} D[i,j]² / |group(g[i])|   # within-group
s_A   = s_T - s_W                                     # between-group
F     = ( s_A / (n_groups-1) ) / ( s_W / (n-n_groups) )
```

`p` is the fraction of permuted `F` values ≥ the observed `F`, with the
standard +1 correction:

```
p = (count{ F_perm ≥ F_obs } + 1) / (n_perm + 1)
```

`n_perm == 0` yields `p = 0.0` ("just to have a deterministic value").

## Call graph

```
skbb::permanova(n, mat, grouping, n_perm, seed, &fstat, &pvalue)   permanova.cpp:351/363
└── permanova_T<TFloat>                                            permanova.cpp:321
    ├── allocate permutted_fstats[n_perm+1]
    ├── permanova_all_T                                            permanova.cpp:277
    │   ├── n_groups = max(grouping)+1 ; group_sizes[n_groups]
    │   ├── permanova_perm_fp_sW_T   -> fills permutted_fstats with s_W
    │   ├── s_T = sum_upper_square(mat)/n                           permanova.cpp:256
    │   └── convert s_W -> F in place
    └── p from permutted_fstats[1..n_perm] vs permutted_fstats[0]
```

`permutted_fstats[0]` is the **unpermuted** case; indices `1..n_perm` are the
permutations. The `s_W` array and the `F` array are the same buffer, rewritten
in place.

## Chunked permutation loop (`permanova_perm_fp_sW_T`, permanova.cpp:68)

```
PERM_CHUNK = <selected backend>::pmn_get_max_parallelism()
step_perms = min(n_perm+1, PERM_CHUNK)
permutted_groupings : uint32[step_perms][n]   # one grouping row per slot
randomGenerators    : RandomGeneratorArray(step_perms, seed)

all rows initialized to the original grouping (once)
for tp in 0, step_perms, 2*step_perms, … < n_perm+1:
    max_p = min(tp+step_perms, n_perm+1)
    parallel for p in [tp, max_p):
        if p != 0: portable_shuffle(row[p-tp], n, randomGenerators[p-tp])
    <backend>::pmn_f_stat_sW(n, mat, max_p-tp, rows, inv_group_sizes, &sW[tp])
```

Three things follow from this shape and matter:

* Chunking exists for **cache locality** (and for the GPU, to bound the device
  buffer). The chunk size comes from the backend, so it is a function of
  `OMP_NUM_THREADS` on CPU and of SM count on GPU.
* Buffer rows are **re-shuffled from their previous permuted state**, not reset
  to the original grouping each chunk. Still a valid permutation, but it means
  permutation `p` depends on every chunk before it.
* Generator `p-tp` is reused across chunks, advancing its state each time.
  Reproducibility therefore requires an identical `step_perms`. See
  [`determinism-and-numerics.md`](determinism-and-numerics.md).

`inv_group_sizes[i] = 1/group_sizes[i]` is precomputed once so the kernel never
divides.

## `pmn_get_max_parallelism()` (permanova_dyn_impl.hpp:51)

| Build | Value | Rationale |
|---|---|---|
| CPU + OpenMP | `2 * omp_get_max_threads() * 16` | 2× threads to amortize spawn, ×16 for the NBLOCK unroll |
| CPU + `SKBB_WASM` | `32` | = 2·1·16, i.e. what native single-threaded gives, so seeded results match |
| CPU, neither | `#error` | deliberate: silently picking a chunk size would silently change seeded p-values |
| CUDA / HIP | `200 * multiProcessorCount` | ~64 blocks/SM to saturate, ×3 for load imbalance |
| OpenACC / OMP-target | `4000` | no portable SM query |

## CPU kernel (permanova_dyn_impl.hpp:114-226)

Two nested pieces of blocking:

* **`NBLOCK = 16`** — 16 *permutations* are evaluated in one pass over the
  matrix, so `mat` is streamed once per 16 permutations instead of once each.
  This is the "better CPU cache locality" work from PR #5.
* **`TILE = 128`** — the (row, col) iteration over the strict upper triangle is
  tiled 128×128 so `grouping_arr[i][col]` stays hot.

The inner loop reads `val = mat_row[col]` unconditionally ("speculatively read,
we will likely use it at least in one of the ifs") and then adds `val*val` to
whichever of the 16 accumulators have `grouping[i][col] == grouping[i][row]`.
Accumulation is in `double` even for the `float` instantiation — deliberate, to
bound accumulation error; the `TFloat` only controls the matrix storage type.

`n_grouping_dims` is decomposed greedily as 16-blocks, then 8, 4, 2, 1.

> **Bug (verified).** The tail decomposition passes the *same* `grouping_arr`
> base to every sub-block instead of advancing it, so the trailing
> permutations of a partial chunk are computed from the wrong grouping rows.
> Triggers when `(n_perm+1) mod 16 ∉ {0,1,2,4,8}`. See
> [`known-issues.md`](known-issues.md#permanova-nblock-tail-uses-the-wrong-grouping-rows).

## GPU kernels

**CUDA/HIP** (`pmn_f_stat_sW_cuda_one`, line 231): one block per permutation,
128 threads, `TILE = blockDim.x`. Each block stages `grouping[row]` for the tile
into `__shared__ row_grouping[128]` (coalesced read; scattered reads then hit
shared memory), accumulates per-thread partials in `double`, then does a 3-stage
tree reduction (128→32→8→1) through `__shared__ double s_W_arr[128]`. Both the
128-thread block size and the two `__shared__` arrays are hard-coded to 128 —
they must be changed together.

**OpenACC / OpenMP-target** (line 330): a much simpler `gang` over permutations
with a `vector reduction(+:s_W)` over columns; no tiling. This is the older,
"experimental" path.

Note the arithmetic differs slightly between paths: the CPU kernel multiplies by
`inv_group_sizes[group_idx]` once per *row* (after accumulating that row's
contributions), while the GPU kernels multiply per *element*. Mathematically
identical, different rounding.

## Complexity and memory

* Time: `O(n² · n_perm / P)` — the matrix is scanned once per permutation.
* Extra memory: `n × step_perms` uint32 for the grouping rows,
  `n_perm+1` TFloat for the statistics, `n_groups` TFloat for the inverses.
  The distance matrix itself is not copied on CPU; on CUDA/HIP a full `n²`
  device copy is made.

## Caller contract (enforced nowhere — validate upstream)

* `grouping` values must be **dense and 0-based**: the code takes
  `n_groups = max(grouping)+1` and indexes `group_sizes` by label. A sparse
  labelling like `{0, 7}` gives empty groups → `inv_group_sizes = 1/0 = inf`.
* `n_groups >= 2`, else `1/(n_groups-1)` divides by zero.
* `n_groups < n`, else `1/(n-n_groups)` divides by zero.
* `mat` must be symmetric with a zero diagonal; only the strict upper triangle
  is read.
* `mat` is `n*n` contiguous row-major; there is no stride parameter.
