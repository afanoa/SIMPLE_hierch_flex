---
name: ole-recipe4-tree-result-2026-10-28
description: "full OLE tree on recipe4 (3+3, crops, embed fix, maps+rho0.03, site_layer=yes): 5.5 h, 18 leaves, ARI 0.15/0.18 to old trees, no leaf recovered, 619 GB without cleanup"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-28T10:02:22.760Z
---

OLE recipe4 tree (seed 1, job 54457724, /scratch/cluster_scratch/atifao/recipe4/trees/ole_recipe4_s1), 2026-10-28 00:23 -> 05:52 (5.5 h).
- 18 leaves (3 under 2,000 rows) vs 31 / 26 in the two older OLE trees (6_, 7_; 15 / 13 small).
- ARI to the old trees 0.145 / 0.178; old trees to each other 0.124. Zero old leaves recovered at Jaccard >= 0.5 (mean best 0.10).
- 619 GB on disk: recipe4 has crops but no cleanup patch; recipe7 adds cleanup (qtoy2: 14 GB removed, 2.9 GB left).

**Why:** confirms the reproducibility problem end-to-end: the new recipe agrees with old runs no better than they agree with each other.
**How to apply:** end-to-end repro test = 1se K switch seed pair (ole-repro session); see [[integration-plan-2026-10-28-builds]], [[split-test-accepts-quality-split-on-null-2026-10-27]].
