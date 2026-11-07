---
name: flex-pca-halfset-fits-lose-fsc-regularization
description: single-halfset flex_pca fits silently run the column solve UNREGULARIZED; the local-eo-split fix regressed out of the tree and was RE-APPLIED 2026-08-14 as cov_needs_local_eo/cov_fit_eo — uncompiled
metadata: 
  node_type: memory
  type: project
  originSessionId: 991559ff-d34f-451c-ad5a-ac1ba95b88c7
  modified: 2026-08-10T16:16:33.003Z
---

`regularize_and_merge_columns` drives a per-shell generalized-FSC Wiener ridge over the covariance
columns — a real, iterated, per-shell prior (`Psh`/`Rsh`, `shrink = <H>/(<H>+R)`). It works in a
normal run. It **collapsed to nothing in every single-halfset fit**: `SIMPLE_COV_BAG=1`,
`heldout=yes`, and the autolp sweep.

Mechanism: inside one halfset every particle carries the same project `eo`, so `Hcol_o` (or `_e`)
stays empty, the `de <= COV_DENSITY_FLOOR .or. do_ <= COV_DENSITY_FLOOR` guard skips EVERY voxel,
`den_e`/`den_o` are 0, the FSC is 0 at every shell, `snr` floors, and the ridge never leaves its
`COV_DENSITY_FLOOR` seed. Signature in the log — all three together:

    fsc= 0.0000   |col|=  0.00E+00   prior=  1.000E+06   shrink= 1.0000

`prior=1e6` is exactly `1/COV_DENSITY_FLOOR` and `shrink=1.0` means nothing is suppressed. The run
still reports success and still produces a basis, so nothing flags it.

Measured on IgG 6k, box_crop=64 lp=12, everything else identical, 2026-08-10:
- no bagging: `fsc` rises 0.36 -> 0.91 by shell 9, `shrink` tracks it (0.36-0.91)
- `SIMPLE_COV_BAG=1`: every shell 0.0000 / 1.0000
- NOT a distributed bug — nparts=5 at 6k reproduces the healthy profile exactly. Bagging is the
  discriminator.

**Fix — regressed out during the Aug-14 thread-over-particles restructure, RE-APPLIED 2026-08-14**
(not yet compiled or run; the user builds). Now two module helpers in `simple_flex_pca_columns.f90`:
`cov_needs_local_eo(build,pinds,nptcls,label)` detects the single-halfset list and logs per stage,
`cov_fit_eo(build,iptcl,local)` returns the halfset. Wired at ALL THREE sites that split on eo —
`cov_stage_subsample`, `accumulate_covariance_columns`, `probe_subspace_iteration` — so no raw
`get_eo` remains in the fitting path. NOTE `cov_fit_eo` is deliberately NOT declared `pure`:
`build%spproj_field` is a `class(oris), pointer`, and a type-bound call through it inside a pure
function was not worth the unverifiable build risk.

The rule, unchanged: `accumulate_covariance_columns` detects a single-halfset particle list and splits
it internally with `eo = merge(0,1,mod(iptcl,4) < 2)`. mod 4, not mod 2, BECAUSE project eo already
alternates by particle index — a mod-2 rule reproduces the project split and leaves the list just as
degenerate. Pure function of `iptcl`, so workers agree without communicating. Normal runs are
untouched (both halves present -> project eo as before). After the fix, bagged runs show
fsc 0.35 -> 0.90 and shrink 0.28-0.90.

**Why this matters beyond the bug:** the entire 2026-08-10 lp sweep in
[[flex-pca-lp-is-inert-at-defaults]] was measured on UNREGULARIZED half-fits, so the lp gain there is
a brick-wall band limit substituting for the missing per-shell prior. Any CV that scores by fitting
disjoint halfsets has this flaw built in — it selects a setting good for an unregularized fit and
then applies it to a regularized production fit. Re-measure before trusting those numbers.

Note the SAME root cause fails differently in the probe M-step: there the FSC is computed against
nothing and the Wiener filter KILLS the basis, so the run dies loudly at the "embedding collapsed"
guard. Column merge = silent under-regularization; probe = loud over-regularization. Documented in
`doc/for_developers/flex_pca_method.md` section 21.

**How to apply:** grep any flex_pca log for `shrink= 1.0000` at every shell — that is the tell that
the column solve is unregularized, whatever else the run reports.
