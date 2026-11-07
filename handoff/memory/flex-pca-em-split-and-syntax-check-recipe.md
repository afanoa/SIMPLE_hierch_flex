---
name: flex-pca-em-split-and-syntax-check-recipe
description: "flex layout after the 2026-09-09 consolidation (22 files; em parent + 9 submodules, states, planes, ICA deleted) and the local no-build validation recipe for flex edits"
metadata: 
  node_type: memory
  type: project
  originSessionId: acaf3fb4-96b2-4601-809e-3c0e25ef0a24
  modified: 2026-09-09T00:00:00.000Z
---

On 2026-08-31 `simple_flex_pca_columns.f90` was split into the procedural parent
`simple_flex_pca_em.f90` plus sibling submodules. On 2026-09-09 (branch flex-pca-resident,
uncommitted at the time of writing) the directory was consolidated from 32 to 22 files:
- `simple_flex_pca_em_{state,env}` -> `_em_fit`; `_em_solve` -> `_em_basis`; `_em_polar`'s 11
  routines -> the module `simple_flex_pca_polar` (its `use simple_flex_gpu` was dead);
  the module `simple_flex_pca_crossfsc` (record types + I/O) -> the em parent.
- `simple_flex_pca_plane_cache` -> `simple_flex_pca_planes` (and `prep_imgs4projected_model`
  moved from latent_ops into planes to break the cycle).
- `gmm + targets + weights` -> `simple_flex_pca_states`.
- `cov_env_int` / `cov_env_int_pub` now live in `simple_flex_pca_util`; the em parent re-exports nothing.
- `simple_flex_pca_ica` DELETED (write-only readout + an env-gated STATE_ICA arm; nothing read it).
Remaining em submodules: fit, basis, crossfsc, mean, embed, estep, mstep, iter, pairmerge.
Public API, CLI flags, env vars, and on-disk filenames unchanged. The same day every flex file had
its notebook-style comments cut (3076 -> 2150 comment lines) with a comment-stripped code-identity check.

**Why:** the flex chain was fully serial and any low-level edit recompiled ~16k lines; the file
count and the lab-notebook comments had grown past what the rest of the repo tolerates.

**How to apply:** validate flex edits without a project build:
`/Users/atifao/.local/bin/devtools-gfortran -fsyntax-only -cpp -ffree-form -fimplicit-none
-ffree-line-length-none -fopenmp -DGNU -DUSE_F08=1 -D__FreeBSD__ -DMACOSX
'-D__FILENAME__="chk.f90"' -I src/inc -I <staged_mods> -I build/modules -J <staged_mods>`.
`-fsyntax-only` DOES write .mod/.smod, so stage fresh flex mods in topo order (util, rounds,
deconv, umap, latent_ops, planes, polar, gpu(-DUSE_FLEX_CUDA), em parent, em submodules, rec3D,
merge, states, model). build/modules is from Aug 31: its `parameters` lacks nfits/l_umap/pindfile
and `load_sigma2_groups` has an old signature, so em_fit and model always show exactly those
errors; run a HEAD worktree first as the baseline, or sed those members in a scratch copy to
unmask anything behind them. New/removed files need a cmake re-run in every build dir (glob lacks
CONFIGURE_DEPENDS). gpu.f90 needs `-DUSE_FLEX_CUDA` to syntax-check.

To go one better and actually RUN a new leaf module's self-test without a project build, compile
just it plus a driver and link against the already-built archive: `build/src/libSIMPLE3.0.0.a` with
`-I build/modules`. The library list to append is in
`build/production/CMakeFiles/simple_test_<id>.dir/link.txt`. Only valid while build/ is newer than
the code the new module depends on. Related: [[flex-pca-lives-as-uncommitted-wip]], [[flex-cuda-build-option]].
