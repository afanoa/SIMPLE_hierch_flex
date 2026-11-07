---
name: per-window-refine3d-shared-reference-inflates-reproducibility
description: TRAP: refining a particle subset with refine3D (both halves against the subset's own map) makes ANY between-subset map difference reproduce at ~0.9 in the coarse bands; random windows did it too; use independent half references
metadata:
  type: feedback
---

RNase P pose test (2026-10-14, `rnasep_coupled/tests/pose`, `tests/pose2`): the two extreme synchronised windows
refined separately with refine3D (refine=neigh, 5 iterations, from their own maps) raised the even/odd reproducibility
of their difference from 0.32/0.39 to 0.94/0.85 at 117-30 / 30-20 A. The null control (random-within-class windows,
no real difference, unrefined 0.15/0.12) went to 0.90/0.87 under the identical protocol.

**Why:** SIMPLE refine3D aligns both half-sets against the same reference (the subset's full map, matching low-pass
~10 A), so the reference's irreproducible noise is imprinted into both halves; the even/odd between-group metric then
reports the reference, not the data. Rigid drift (0.2-0.5 deg, 3-5 A per subset) adds a reproducible gradient term in
the fine bands on top; always re-register refined maps before differencing.

**How to apply:** never read a between-state / between-window gain off maps whose halves were refined against a shared
reference. Split each group into disjoint halves, reconstruct and refine each half from its OWN map, and correlate the
A-run difference with the B-run difference (gold standard), or compare against the random-window null refined the same
way. Related: [[view-imbalance-inflates-between-state-metric]], [[rnasep-2d-fine-3d-coarse-complementary-2026-10-14]].
