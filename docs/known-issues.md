# Known issues and latent hazards

Findings from a full read of the tree. Each entry says how it was established —
**verified** (reproduced here), **read** (from the source, not executed), or
**latent** (works today, will break under a plausible change).

Nothing in this file has been fixed. Do not "clean up" any of it as a drive-by
inside an unrelated change — raise it, get agreement, fix it in its own commit
with a test.

---

## PERMANOVA: NBLOCK tail uses the wrong grouping rows

**Verified.** `src/distance/permanova_dyn_impl.hpp:186-224`.

`pmn_f_stat_sW_cpu` builds one array of 16 grouping-row pointers per 16-permutation
block, then decomposes a partial block greedily as 8+4+2+1. The sub-block calls
advance the *output* offset (`gblock2`) but always pass the **same**
`grouping_arr` base, so every sub-block after the first re-reads grouping rows
`gblock+0…` while writing to `group_sWs[gblock2…]`:

```cpp
uint32_t gblock2=gblock;
if ((gblock2+(8-1))<n_grouping_dims) {
  pmn_f_stat_sW_block<TFloat,8>(n_dims,mat, grouping_arr, …, group_sWs+gblock2);
  gblock2+=8;
}
if ((gblock2+(4-1))<n_grouping_dims) {
  pmn_f_stat_sW_block<TFloat,4>(n_dims,mat, grouping_arr, …, group_sWs+gblock2);
  //                                        ^^^^^^^^^^^^ still rows gblock+0..3
```

Effect: the affected permutations get a duplicate of an earlier permutation's
`s_W` instead of their own. The p-value is then computed from a permutation
distribution containing duplicates. When the *first* chunk is the partial one
(i.e. `n_perm+1 < PERM_CHUNK`), the duplicated row can be permutation 0 — the
**unpermuted** statistic — which biases `p` upward.

Trigger condition, exactly: **`(n_perm+1) mod 16 ∉ {0, 1, 2, 4, 8}`**, i.e. 11
of every 16 values of `n_perm`. It does not depend on the thread count —
`PERM_CHUNK = 2·OMP_NUM_THREADS·16` is always a multiple of 16, so every chunk's
size is congruent to `n_perm+1` mod 16.

Every `n_perm` used anywhere in the test suites lands on a safe residue, which is
why this was never caught: 999 → 8, 99 → 4, 199 → 8, 499 → 4. `n_perm = 9`
(→ 10) and `n_perm = 1000` (→ 9) both trigger it.

Reproduction: a 2×2 distance matrix makes the expected `s_W` per grouping row
trivially known (1.0 if the two samples share a group, else 0.0). Pass 10 rows —
eight `{0,0}` then two `{0,1}` — and rows 8 and 9 come back as 1.0:

```
row 7  grouping {0,0}  s_W got 1.0  want 1.0  ok
row 8  grouping {0,1}  s_W got 1.0  want 0.0  WRONG
row 9  grouping {0,1}  s_W got 1.0  want 0.0  WRONG
```

Full reproducer in the tracking issue.

Likely fix: pass `grouping_arr + (gblock2 - gblock)`. Any fix changes seeded
p-values for the affected `n_perm` values, so the WASM expected headers must be
regenerated — see [`determinism-and-numerics.md`](determinism-and-numerics.md#things-that-silently-change-seeded-results).

Related, same function: `grouping_arr[i] = groupings + (gblock+i)*n_dims` is
computed for all 16 `i` even when `gblock+i >= n_grouping_dims`. The out-of-range
pointers are never dereferenced, but forming them is UB.

## PERMANOVA GPU: `permutted_sWs` device buffer is one element short

**Read** (no GPU available here to execute).
`src/distance/permanova.cpp:135, 146, 232, 241`.

The host buffer holds `n_perm+1` statistics (allocated in `permanova_all_T` as
`permutted_fstats[n_perm+1]`, and the loop writes indices `0..n_perm`), but the
device allocation and copy-back both use `n_perm`:

```cpp
skbb_acc_nv::acc_create_buf(permutted_sWs, &permutted_sWs_device, n_perm);
…
skbb_acc_nv::acc_copyout_buf(permutted_sWs, permutted_sWs_device, n_perm);
```

Under CUDA/HIP that is a one-element device heap overflow on write, plus the last
permutation's `s_W` is never copied back (the host slot keeps whatever
`new TFloat[]` left there). Under OpenACC/OpenMP-target the last element is
outside the mapping.

## `api_tests` exit code masks earlier failures

**Verified by reading.** `api_tests/test_distance.c:124`,
`api_tests/test_ordination.c:436`.

Both `main()`s end with `return failed ? 1 : 0;`, but `failed` is reset to `0` at
the top of every test function. A failure in an earlier function is therefore
invisible to the exit code if a later one passes. `global_failed` is maintained
for the printed summary but never used for the return value. CI's `make test`
consequently cannot fail on those cases.

## Test asserts the wrong variable after an fp32 call

**Verified by reading.** `src/tests/test_permanova.cpp:54-59` (and the same
pattern at `api_tests/test_distance.c:60-65`).

```cpp
skbb::permanova(n_samples, matrix_fp32, grouping_equal, 999,
                rand_seed_invalid, stat_fp32, pvalue_fp32);
ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);      // fp64, from the PREVIOUS call
ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);
```

The fp32 result of that call is never checked.

## Stale comments that contradict the code

**Verified by reading.** `src/tests/wasm/generate_permanova_expected.cpp:15-25`
says the generated header "is checked into the repo so the WASM test is
hermetic" and tells you to regenerate with
`make -C src wasm_regen_permanova_expected`. Both are false: the headers are
`.gitignore`d (that was the maintainer's PR #12 request) and no such make target
exists. `src/wasm/emscripten_build.mk:120-125` states the opposite, correctly.

This is the failure mode the comment rules in `CLAUDE.md` are about: the comment
described the state of one pull request rather than the code.

## Unguarded `PREFIX` in `install` / `clean_install`

**Read.** `src/Makefile:394-403`. `PREFIX` falls back to `CONDA_PREFIX`; if both
are unset the recipe runs `mkdir -p /lib` and copies there. `install_inmem` in
`src/inmem_build.mk:70-74` has the guard the others lack.

## `ordination.h` documents the output shape ambiguously

**Read.** `src/extern/ordination.h:104` describes `samples` as
"Matrix of size (n_eighs x n_dims)" while the adjacent prose says "row-indexed by
the sample id". The code produces `n_dims` rows of `n_eighs` values. The element
count is identical, so nothing breaks — but a caller indexing
`samples[axis*n_dims + sample]` gets garbage.

## Latent build fragilities

* `make -C src clean` runs `rm -f *.cpp *.hpp *.h *.cu` in `src/`. It is safe
  only because every such file there is generated. It also misses `*.hip`
  (generated by the AMD HIP path).
* `test_permanova.exe` links `libskbb_cpu.a` **without** `$(BLASLIB)`. It works
  because static-archive linking pulls only the objects it needs and PERMANOVA
  touches no BLAS. Any new dependency from `distance/` onto `ordination/` breaks
  that link.
* `SKBB_OBJS`, `WASM_OBJS` and `INMEM_OBJS` are three hand-maintained lists of
  the same translation units. Nothing checks that they agree.
* No top-level `inmem_static` target; it is reachable only as
  `make -C src inmem_static`.
* `src/util/skbb_dl.cpp:31,74` tests `SKBB_GPU_INFO[0] == 'Y'` while
  `skbb_detect_acc.cpp` uses a negative-list convention, so `SKBB_GPU_INFO=y`
  enables one set of messages and not the other.
* CUDA/HIP failures `throw std::runtime_error` through an `extern "C"` boundary;
  the public API has no error channel and nothing catches them.
