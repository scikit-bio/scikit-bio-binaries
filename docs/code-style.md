# Code style and change discipline

This repository is upstream `scikit-bio/scikit-bio-binaries`, maintained by
Igor Sfiligoi (@sfiligoi). We contribute through PRs from a fork. **Conformance
to what is already there beats our taste**, and every line we touch is a line the
maintainer has to review.

## There is no formatter — so churn is never free

No `.clang-format`, no `make format`, no format check in CI. Nothing will
reformat a file back, and nothing enforces a house style. Therefore:

* **Never reformat, re-indent, re-wrap, or re-order anything you are not
  functionally changing.** A whitespace-only hunk is pure reviewer cost.
* Do not fix typos in comments or identifiers you are not otherwise editing.
  The codebase contains `permutted`, `avaialble`, `cpecific`, `trashing`,
  `unction`, `n_eighs` (for eigenvalues), `skbio_bins` in log strings. They are
  harmless and some are load-bearing (`n_eighs` is in the public API).
* Do not "modernize" (`typedef`→`using`, `malloc`→`new`, raw loops→algorithms)
  outside the scope of the task.

## Match the file, not a global rule

Indentation genuinely varies by file and there is no correct answer other than
the surrounding code:

| File | Indent |
|---|---|
| `src/util/skbb_detect_acc.cpp` | 1 space |
| `src/distance/permanova.cpp`, `src/ordination/principal_coordinate_analysis.cpp` | 2 spaces, with tabs in some continuation lines |
| `src/ordination/linalg_backend_*.cpp`, `src/util/portable_shuffle.hpp`, `src/tests/wasm/*` | 4 spaces |
| generated `.cpp`/`.h` | whatever the generator emits — never hand-edit |

Other conventions that are consistent and should be followed:

* BSD 3-Clause header block at the top of every source file. Files derived from
  UniFrac carry both copyright lines (`2016-2025, UniFrac development team` and
  `2025--, scikit-bio development team`); new files carry only the latter.
* `template<class T>`, not `template<typename T>`.
* Type parameter names: `TFloat` (float/double), `TNum` (numeric+bool),
  `TReal`/`TRealIn` in the ordination code.
* Implementation functions carry the `_T` suffix — **this is parsed by the code
  generators**, see [`codegen.md`](codegen.md).
* `snake_case` for functions and variables; `n_dims`, `n_perm`, `n_eighs`,
  `n_groups` for sizes.
* `} else {` on one line; no `using namespace`.
* Public headers are C-compatible: `EXTERN` macro, `#include <stdbool.h>` under
  C, no C++ types in signatures.
* Memory: both `new[]/delete[]` and `malloc/free` are in use, sometimes in the
  same file. Match the function you are editing rather than converting it.

## Comments

The existing comments are short and explain **why**, not what:

```cpp
for (uint32_t col=row+1; col < n_dims; col++) { // diagonal is always zero
// Use full precision for intermediate compute, to minimize accumulation errors
// speculatively read, we will likely use it at least in one of the ifs
```

Write comments that will still be true and useful to someone reading this file in
two years with no knowledge of the change that introduced them.

**Long comments are justified only when they encode a constraint that would
otherwise be rediscovered the hard way.** Good examples already in the tree:

* `util/portable_shuffle.hpp` — why `std::shuffle` is banned.
* `linalg_backend.hpp` — the `svd_no` stride contract.
* `permanova_dyn_impl.hpp:64-70` — why the `#error` exists instead of a default.
* `Makefile:5-9` — why `.NOTPARALLEL` is required.

**Do not narrate the work.** These are the failure mode, and the tree already
contains examples that have gone stale:

* `"Stage 3: …"`, `"Stage 6: …"` — PR-plan numbering in source files.
* `"per reviewer guidance"`, `"(see PR #12 discussion)"` — review-round context.
* `"Exact behaviour of the pre-refactor principal_coordinate_analysis.cpp"` —
  meaningless once nobody remembers the refactor.
* `generate_permanova_expected.cpp:15-25` — describes a policy (headers checked
  in, a `wasm_regen_*` make target) that was **reversed during review**; the
  comment now actively misleads. See
  [`known-issues.md`](known-issues.md#stale-comments-that-contradict-the-code).

That context belongs in the commit message and the PR description, where it is
timestamped and does not rot in place.

## What the maintainer has actually pushed back on

From the review history — treat these as standing requirements:

1. **Duplicated build infrastructure.** On PR #12: *"Most of this file looks
   exactly like the Makefile. Couldn't we consolidate the two, so we do not need
   to maintain two copies? Especially as/when we add additional functionality."*
   → produced the `skbb_cpu_tu` canned recipe. Add a new translation unit in
   **one** place; do not copy a rule between `src/Makefile`,
   `wasm/emscripten_build.mk` and `inmem_build.mk`.
2. **Committing generated artifacts.** On PR #12: *"Why do we check in the
   expected files, if they have been dynamically generated?"* → the WASM
   expected-value headers are now `.gitignore`d and regenerated on demand.
   Generated output does not go in git.
3. **Complexity is tolerated but noticed.** On the GPU/x86 dispatch rules:
   *"A bit complex, but cannot think of a better way."* Prefer the boring
   solution; if a change makes the build or dispatch harder to follow, say so in
   the PR description and explain why the simpler option does not work.

## Commits and PRs

* Stage explicit paths — never `git add -A` / `git add .` / `git commit -a`. The
  tree fills with build artifacts (`*.o`, `*.a`, generated `*.cpp`) that are not
  all gitignored, and they will end up in the commit.
* Confirm before pushing: `git status`, then
  `git diff --name-status <base>...HEAD`.
* Keep the diff to the task. Unrelated fixes — including the entries in
  [`known-issues.md`](known-issues.md) — go in their own commit with their own
  test.
