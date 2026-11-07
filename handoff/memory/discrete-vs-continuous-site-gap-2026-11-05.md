---
name: discrete-vs-continuous-site-gap-2026-11-05
description: "IgM-like fixture: bin-map rank/affine tests cannot certify discrete site states at realistic d' (~1); 1-D path dips are fragile; GAPw (consensus density between cross-half first-factor lobes) separates a1 (fixture .02-.23 at >= 50k vs null >= .37, 4 arms), misses 3-in-a-row a2"
metadata:
  node_type: memory
  type: project
  originSessionId: be3d41f9-61a0-48c1-b522-ed1c316c2df1
  modified: 2026-11-06T09:46:48.404Z
---

2026-11-05, side session of the IgM fixture study. Work dir: batch fixture_igm/discrete/. Write-up: discrete_vs_continuous.md in the
shared notes folder of session local_c1bafd06.

- Oracle (truth-built bins, SIMPLE reconstruct3D box_crop 128): the null's continuum signature (second bin-map factor) has z 13-23
  with exact-angle bins but 8.5 / 1.9 / 0.3 at d' 4 / 2 / 1.1. Real site classification reaches d' ~1.0-1.1, so rank /
  affine-residual tests can only reject two states, never confirm them.
- 1-D consensus dip along the first-factor path: fixture .11-.18 vs null <= .04 on oracle paths. A 10 deg path tilt gives false
  null dips in up to 4/12 directions, and input-free paths differ 12-18 deg between halves. Rejected.
- GAPw: lobes of w = sign(F1e+F1o) sqrt(max(F1e F1o, 0)) at 0.5 of the extreme; B = their convex hull minus the lobes dilated
  10 A; consensus density in B / the lower lobe density.
  - Gates: r1 >= .5 (signal-dominated factor), lobe half-separation >= 10 A, even/odd lobe axes <= 30 deg, |B| >= 10 voxels.
  - Fixture a1: .02-.09 (oracle), .20 / .23 / .32 (input-free at 100k / 50k / 25k). a3 / a4 two-state arms (.65/.35, .70/.30,
    40 deg): .13-.15 (oracle). Null >= .37 on 26 sets (four arms).
  - Threshold .28 = null median - 3 robust sd. A lobe threshold of 0.4 gives false null gaps on weak orderings: use 0.5.
- Needs a cross-validated ordering: each half scored against the other half's site direction, information-stratified bins.
  Merged-reference stage-B bins of the null carry a reproducible noise-locked factor (z1 38.9).
- Misses a2 (three states in a row: the middle state fills the gap). Sensitive to the lobe threshold: only 0.5 has margin.
- 2026-11-06, the xstageb v2 swap battery (14 runs, `discrete/scripts/dc_battery.py`):
  - fixture a1 DISCRETE in both seeds (.08-.10); 0 / 9 null and core runs DISCRETE;
  - the rigid-core patch controls pass G1 (cross-half r .9), and GAPw stops them (1.45-1.54);
  - most arm sites give no call (between-lobe region < 10 voxels).
- 2026-11-06 gate replay (xsbg.py runs; dataset-generic call `discrete/scripts/dc_gates.py`, grid = nearest even box to 5.11 A/px):
  - 0 / 56 DISCRETE: qtoynull 20, olenull 4, RNase P 32 (incl. the 692k closed leaves); v_fix_a1 DISCRETE .09 as validation;
  - RNase P per-site stage-B splits never separate the closed reference rows, so the layer adds no split.
- 2026-11-06 wave 2 (33 runs on the 692k closed lineage's other nodes r13 / r131 / r1311 / r132 / r21 / r211):
  - 0 DISCRETE, and SEPARATES no throughout;
  - lf_r211_s01 is real (r1 .98) but not discrete (GAPw .45);
  - lf_r21_s01-03 split open vs intermediate and are stopped at G2 / G3 (lobes 6-8 A apart).
- TRAP: node stageB/b.simple refers to the full project's stacks, which SIMPLE selection rejects; use the prep's node.simple.

**Why:** the coordinating session needs a discreteness gate, so the tree splits discrete arms but not the null twin's continuous swing.
**How to apply:** GAPw is a prototype verdict on this one fixture. A continuum's GAPw is set by the fan's geometry, so a real node needs
its own null level. Related: [[tree-good-pose-null-fails-halves-and-seeding-2026-10-30]], [[soft-multistate-fixed-pose-verdict-2026-10-30]],
[[tree-continuum-fixture-and-verdict-2026-10-29]].
