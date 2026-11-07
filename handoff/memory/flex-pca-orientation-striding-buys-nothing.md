---
name: flex-pca-orientation-striding-buys-nothing
description: Orientation-stratified subsampling in flex_pca is inside the noise floor on IgG; the proportional mode doubles as a resampling control that puts that floor at 0.05 mean cosine
metadata: 
  node_type: memory
  type: project
  originSessionId: 991559ff-d34f-451c-ad5a-ac1ba95b88c7
  modified: 2026-08-10T16:16:51.551Z
---

`SIMPLE_COV_ORISTRIDE` (added 2026-08-10, default OFF) stratifies `cov_stage_subsample` by
projection-direction cell: 1 = proportional (reproduces the full angular distribution exactly),
2 = equalized (water-filled, flattens preferred orientation).

**It buys nothing on IgG.** Measured, box_crop=64 lp=12, 30k particles, caps rescaled to production
binding ratios, cross-halfset basis cosine:

| mode | mean cos | n>0.9 |
|---|---|---|
| 0 index stride (current) | 0.5799 | 1 |
| 1 proportional (control) | 0.5301 | 0 |
| 2 equalized | 0.5367 | 0 |

**The control is the reusable part.** Mode 1 has an IDENTICAL angular distribution to mode 0
(chi2/dof 0.00 against the full set) but shares only 25% of its particles, so mode0-vs-mode1 is a
pure resampling control. It puts the pipeline's noise floor at **0.050 mean cosine**. Mode 2 sits
0.007 from mode 1 — an order of magnitude inside that floor.

Use that 0.05 floor to judge every other flex_pca A/B scored this way. It already retired one claim:
the bc48-vs-bc64 difference at matched band (0.024) is noise, so "a larger box at matched band is
slightly better" was withdrawn.

**Why there was nothing to fix on IgG:** row order carries no orientation structure there — per-stack
mean direction 0.0278 against 0.0302 for a shuffled control, and every stage's index-strided
subsample reproduced the full angular histogram at chi2/dof < 1. The feature is a guard for projects
ordered by micrograph/tilt/optics group, where a stride CAN alias against orientation. Do not enable
it expecting a gain.
