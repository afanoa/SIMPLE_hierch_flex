---
name: rnasep-goal-is-conformations-not-periphery-resolution
description: "2026-10-03 user decision for RNase P: stop chasing periphery resolution (needs multibody refinement, which they won't pursue); the goal is the conformations, ideally from an automated molecule-in, states-out pipeline"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-03T21:45:01.283Z
---

On 2026-10-03 the user decided the RNase P work should target the conformational states, not peripheral resolution.

Evidence behind the decision:
- Region-1 thirds of S4 reproduce (1−corr 0.19–0.20).
- Their periphery local resolution is no better than a random third: median 6.1–6.8 Å vs random 6.2 Å, consensus 5.7 Å.
- Weighting, robust-frame, joint-pose, solvent-prior and body-frame attempts also gave no periphery gain.

**Why:** the user thinks peripheral resolution is out of reach without multibody refinement, which they consider inelegant and won't pursue.

**How to apply:**
- Measure success by reproducible, distinct conformations: cross-half gap test, selection-bias null.
- Don't measure success by local resolution.
- The user wants it automated: molecule in, states out. The chain is whole-molecule fit → automatic regions (eigengap, port regions_auto to Fortran) → joint region-block fit → reproducible regions only → state count by the gap test (island vs continuum) → refine3D_auto per state → composite.

Related: [[rnasep-multibody-region1-split-is-real]], [[multibody-flex-means-region-pca]], [[flex-pca-state-reproducibility-gap-test]].
