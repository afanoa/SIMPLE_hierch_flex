---
name: rnasep-per-window-poses-gold-standard-2026-10-14
description: RNase P pose test (gold standard, independent half references): refining poses per window along the synchronised 2D coordinate raises the reproducible window difference to 0.72/0.32/0.41 at 30-20/20-12/12-8 A (random null 0.35/0.08/0.05); coarsest band at the random floor; consensus poses are what blur the 2D motion out of 3D
metadata:
  type: project
---

Result (2026-10-14, `rnasep_coupled/tests/pose2`, `tests/pose3`, notes section 19.3, `viz_pose_gold.png`): balanced
synchronised windows 1 and 8 split by member parity into halves A/B, each half reconstructed and refined from its own
map (refine3D refine=neigh, 5 iterations), refined maps re-registered on their starts; A-run difference vs B-run
difference per band 117-30 / 30-20 / 20-12 / 12-8 A:
- synchronised windows: unrefined 0.25 / 0.31 / 0.06 / 0.15 -> refined 0.55 / 0.72 / 0.32 / 0.41
- random-within-class windows (same protocol): 0.17 / 0.19 / 0.02 / 0.01 -> 0.59 / 0.35 / 0.08 / 0.05
The refined difference is the unrefined one made clean (same arm dipole, 17 A apart, lobes 4+3 %, still distributed).

**Why:** consensus poses mispose particles whose structure differs along the 2D coordinate; the misposing blurs the
window maps into each other. Window-conditional poses recover the difference below 30 A. The coarsest band is a
random-composition effect that refinement amplifies for any window set.

**How to apply:** "moving the 2D motion into 3D" needs states/windows with their own references and poses, not a
flex latent on consensus poses. Read every refined window set against the refined random null and only with
independent half references ([[per-window-refine3d-shared-reference-inflates-reproducibility]]). Builds on
[[rnasep-2d-fine-3d-coarse-complementary-2026-10-14]] and [[rnasep-subclass-window-trajectory-2026-10-14]].
