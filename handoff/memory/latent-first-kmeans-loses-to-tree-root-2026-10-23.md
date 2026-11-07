---
name: latent-first-kmeans-loses-to-tree-root-2026-10-23
description: "Latent-first node (global k-means on flex / fflex raw MAP z -> K maps -> E-step) never beats hetero_tree's own root split: CAHRA 0.18-0.24 vs 0.463 (ceiling 0.531), cdisc A 0.10-0.17 vs 0.413; the latents have no stable clusters (ARI <= 0.16); NMI-based stability is biased to large K"
metadata:
  node_type: memory
  type: project
  originSessionId: a4258345-b249-48ff-a0d3-d33e99ae33ac
  modified: 2026-10-24T02:03:11.555Z
---

2026-10-23. A user-approved prototype that the user extended into a development loop. Scripts: batch2
`rnasep_klem/fast_tree/lfn/` (lfn_cluster.py, lfn_score.py, lfn_screen.py; jobs fast_lfn_v1 / fast_lfn_screen_v1).
Roots used: CAHRA tree_m8c36_c11_s1, toys fulltoy_p27cd1 (cdisc) and p27igm.

- Iteration 1, as specified: raw z, K by NMI restart / half stability, reconstruct3D nstates=K, then one refine3D pass with
  stage-A keys as the E-step. The current code has no fixed-pose state scorer: the Oct 5 score_states lived in the softem branch
  only. After the E-step: CAHRA open NMI 0.128, cdisc A_domain 0.124, IgM combo 0.246.
- TRAP: pairwise NMI between partitions rises with K even for chance partitions, so "most stable K" picked K=10 (the top of the
  range) on every set. Use a chance-corrected measure (ARI) for K stability.
- Iteration 2, Python only: {fflex, flex} x {raw z, z / site contrast, z / contrast with c^2 q_jj weights}, with ARI-based K.
  Before the E-step:
  - CAHRA: 0.178, and 0.236 even at the truth-picked best K, against the tree root's 0.463.
  - cdisc: A_domain 0.10-0.15 against 0.413.
  - IgM: 0.165-0.18 against 0.135. That is K-confounded: 8-10 clusters vs a 2-state root split.
- The ARI stability of every latent is <= 0.16 at every K: there are no stable clusters. The state signal is low-dimensional
  (the CAHRA open state sits mainly on z4), so global 10-D k-means is driven by other directions.
- The E-step improves the clusters it is given but cannot create the right partition (as on Oct 5).
- The loop was stopped after two iterations without gain. No port spec.

**How to apply:** don't re-propose global-latent clustering as a replacement for the tree's site / seeds steps. The one idea the
evidence still supports is a per-direction mixture chosen by cross-half reproducibility, which generalises site_rank's direction
scan. Related: [[flex-placement-not-latent-loses-small-states]], [[flex-latent-density-basins-no-power-2026-10-22]].
