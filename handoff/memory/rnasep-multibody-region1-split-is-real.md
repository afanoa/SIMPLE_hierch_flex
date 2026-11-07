---
name: rnasep-multibody-region1-split-is-real
description: "2026-10-01 RNase P S4: multibody region-1 slices refined with refine3D_auto differ 0.20 inside region 1 (no reference bias, noise floor 0.0003); marginal pose polish under the multibody model changes nothing (R 0.073 vs jitter 0.075)"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-01T07:58:58.499Z
---

S4 (238k, hires poses), multibody lp=8 fit `rnasep_S4couple/multibody_lp8` (G10cov_labels_B3, 6 cols/region).
Hard equal-count thirds along region 1's leading column (z1 = coordinates col 6), refine3D_auto per slice
(`~/rnasep_mbauto.sbatch`, T=1/3; T=0/9 = disjoint random thirds; SAMESTART=1 starts from the S4 consensus).

Region-masked 1-corr at 15 A (scorer `~/rnasep_r1score.sbatch`), region 1 / outside region 1:
- multibody flex maps slice 1 vs 3: 0.074 / 0.0015
- refine3D_auto, own-start maps: 0.202 / 0.005 (3.50 / 3.55 A, 79k each)
- refine3D_auto, SAME start: 0.202 -> no reference bias
- random third vs random third: 0.0003 -> noise floor; all at 3.44-3.55 A
- for scale: S4 whole-molecule substates 0.07-0.08 in region 1 (leaking into region 3); S1 vs S4 0.127
The user sees a domain appearing/disappearing in region 1. The flex delivery maps understate it ~3x.

Pose polish (step 1 of coupled pose/heterogeneity): combined build `~/SIMPLE_regpose` (local worktree
`~/SIMPLE_regpose`, branch flex-region-pose = flex-region-pca-in2 + the uncommitted flex-pose-marginal polish
patch), SIMPLE_COV_POSE_POLISH=1 DEG=3 PASSES=3 STEP=soft DEFLATE=1 SHIFT=0. Polish moved median 0.53 deg.
Refit with region2 build: polished R 0.0726, jitter control J (same angles, random axes, `~/tools/jitter_poses`)
0.0746, baseline 0.0736 -> null. Hires S4 poses are already good; pose is not the limiter here.

**How to apply:** don't spend more on pose coupling for S4; the lever is localisation (region blocks) plus
hard selection + refine3D_auto. Related: [[flex-region-pca-branch-state]], [[multibody-flex-means-region-pca]].

**CORRECTION 2026-10-03:** "no ref bias / noise floor 0.0003" did not control for SELECTION bias. Thirds along a random spectrum-matched direction give 0.132, so about 2/3 of the 0.20 is bias. The real remainder is a coherent arm change. See [[rnasep-thirds-selection-bias-null]].
