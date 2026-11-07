---
name: joint-pose-heterogeneity-toy-verdict
description: "2026-10-02 toyJ: heterogeneity biases consensus poses (toy3 4.4 deg, toy2c 7.7 deg, moving-region maps 0.75 vs 0.99) but a SOLVAR-type z-marginal discrete pose search only removes ~15% of it; oracle search reaches 1.44 deg, so model fidelity is the limit"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-02T14:49:44.563Z
---

Prototype: numpy EM on batch2 `~/toyJ/` (toyj_core.py conventions resolved by toyj_conv.py; toyj_em.py =
z-marginal score, discrete local pose search, posterior-mean pose anchor, shared-memory M-step, flags
--u_support/--mixture/--oracle; toyj_maps.py = oracle-binned map quality; pose_bins.py = error vs |theta|).
Data `/mnt/beegfs/elmlund/afan/rnasep_regions/toyJ/{s0.02 (toy2c 125 states), t3s0.02 (toy3 31 levels)}`,
40k, SNR 0.02, no astigmatism, est poses from refine3D_auto (`est_ptcl3D.txt`).

Findings:
- refine3D_auto poses are biased by conformation (toy3: 1.4 deg near the mean -> 7.3 deg at |theta|=15; toy2c
  7.7 deg median, no global frame offset). Latent R^2 barely drops (0.95 -> 0.91), but region maps at the swing
  extremes fall 0.99 -> 0.75 (8 A). Pose bias is invisible to latent metrics; score maps.
- Joint search with a whole-molecule basis: null (basis absorbs rigid rotation = gauge degeneracy).
- Basis restricted to the moving region (multibody support): extremes 7.25 -> 6.08 deg, maps 0.75 -> 0.78.
- Rank 12 tight support / soft state x pose mixture (K=6,10): ~3.8 deg median, states collapse (pops 4, 13, 41).
- ORACLE (search scored against each particle's true volume): 4.4 -> 1.44 deg everywhere in 2 iterations.
  The discrete search works; the linear / mixture conformational model is not faithful enough for a +-15 deg swing.
- toy2c has ~no static core (8.5% of dense voxels), so the true frame is unidentifiable there; use toy3.

**How to apply:** the joint idea only pays with a model that represents the moving body exactly, i.e. per-body
rigid poses (multi-body refinement) searched by the same discrete machinery, core pinning the global pose.
Related: [[toy3-body-frame-test-system]], [[solvar-marginalizes-z-during-fitting]], [[rnasep-multibody-region1-split-is-real]].
