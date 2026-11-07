---
name: hetero-tree-toy-root-flex-nondeterminism-and-oneflex-2026-10-22
description: "child splits are start-dependent (r1 child labels NMI ~0.5 same parent, 0.1-0.3 across); ctrue toy: 1 of 3 identical-setting roots diverged in flex_pca (fflex latents |corr| ~0.5-0.8, seeds NMI 0.71, 2-leaf tree); oneflex=yes kept the reference's 5-leaf tree, closed 0.48/0.83 vs 0.52/0.80, children -41 %"
metadata:
  node_type: memory
  type: project
  originSessionId: a4258345-b249-48ff-a0d3-d33e99ae33ac
  modified: 2026-10-23T03:11:08.435Z
---

2026-10-22 night, fast-tree session. Build fastc1 = 51c8cdac4 + threshold cache + oneflex. Runs: `rnasep_klem/fast_tree/toy_ctrue_{ctl,of}`
and `toycmp`; the reference is the port's `stree_port/fulltoy_ctrue` (build 49). The root runs identical settings in all three runs,
since oneflex acts only at depth >= 1.

**Root reproducibility:**
- Two roots agree: oneflex run vs reference, seeds NMI 0.978, fflex latents |corr| 0.94 / 0.91 / 0.85.
- The control diverged: 3 poles instead of 2, seeds NMI 0.71, a 2-leaf tree that misses the domain-A split.
- Between any two runs, flex labels NMI is 0.48-0.68 and fflex labels NMI 0.12-0.25.
- The divergence starts in flex_pca's latents, not in site_rank: site_rank is deterministic given the mask (NMI 1.000).

**oneflex** (children skip the mean, whole flex and cheap steps; parent mask and maps):
- Same topology and 5 leaves as the reference, and the same delivered combinations.
- Closed leaf 0.481 / 0.828 vs 0.520 / 0.795; domain-A-off leaf 0.44 vs 0.28 A; junk identical.
- Child nodes 116 vs 198 min summed; 6-11 min saved per child on the same build.

**Seed test and replicate (2026-10-23):**
- The same SIMPLE_SEED reproduces the root (seeds NMI 0.977), so the root is start-dependent with two attractors.
- 2-pole root (13.2k + 161k): the big child splits and domain A is recovered (5 leaves).
- 3-pole root: the big child merges its endpoints into one pole and stops (2-3 leaves, no domain-A split).
- This holds with and without oneflex, so oneflex changes no outcome. The root's start decides tree quality, which argues for a
  multi-start root.

**Real data, the r1 children** (2026-10-23, the port's m8c16 / m8c17, `fast_tree/wval`):
- The r1 child (54k) splits every time, but its stage-B labels NMI is only 0.49-0.57 between two hard runs from the SAME parent,
  and 0.11-0.29 across parents that agree on 99 % of members.
- Same pathology as the toy: the seeds are one pole holding nearly all particles plus a tiny chance pole (300-870).
- Margin weighting can't fix it (inconclusive at 0.29 vs 0.11 / 0.20). The candidate fix is a pole guard (a non-empty minority
  hard class) plus a data-driven fallback.

**All 10 flex directions at the RNase P r1 child** (`fast_tree/wval/sdirs`):
- 15-19 of 20 endpoint EMs COLLAPSE (hard class ~98 % of the node) yet pass the halves test. The tree's tiny second pole is their
  leftover. The rest end at <= 547 particles.
- So it is a no-signal node, the mirror of the toy's EMPTY classes.
- Pole guard (revised): a pole's hard class must be < 1/2 of the node AND pass the cross-half test on its hard-class maps. Otherwise
  fall back to the next site, then an fflex 2-state start, then a leaf. A per-half size floor of 6 is too weak.

**Why:** the user requires whole landscapes to reproduce across independent runs ([[landscape-must-be-reproducible-and-splits-data-derived]]).

**How to apply:** judge tree changes on at least two independent runs; check fflex latents first when trees differ.
Related: [[hetero-tree-schedule-cuts-fail-unconverged-drift-2026-10-22]].
