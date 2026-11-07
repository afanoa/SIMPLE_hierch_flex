---
name: klem-leaf-recipe-score-coordinate-2026-10-20
description: "leaf recipe for continuous motion at realistic SNR (simswing30): K=8 hard refine3D_states only as basis -> score every particle vs the K maps, PCA of score vectors, rebin/rescore 2-3x -> 10 equal-width leaves on the motion PC; oracle 108k corr 0.735 (labels 0.47, true-map ceiling 0.753), ends 5.1/25.3 A; realistic 75 %-pure family ends 5.0/24.5; boards, K>8, recursion fail; rnasep_klem/lrec30/RECIPE.txt"
metadata:
  node_type: memory
  type: project
  originSessionId: 9e800932-bfcf-4aa9-9843-660dc7ddf3af
  modified: 2026-10-21T05:23:06.838Z
---

2026-10-20/21 (task: minimal leaf-discovery recipe for continuous motion at realistic SNR; benchmark simswing30, swing 0..30 A
of conformers 200-230; see [[klem-bottom-up-organization-2026-10-20]]). Full table: rnasep_klem/lrec30/RECIPE.txt.
- Mechanism: refine3D_states labels by argmin of best-pose euclid dist; HARD labels are the bottleneck (K=4/8/16/20 eta2
  0.22/0.25/0.20/0.17; even TRUTH-seeded K=8 decays to 0.47). Its exec dir keeps the last iteration's tables
  dist_part_neigh_<part>.dat (header int64[4] nrefs nptcls nnz nchunks; chunks at the END: chunk start = size - 28 nnz -
  4 nchunks; 24-byte candidates iref/inpl/dist/x/y/has_sh) -> per-particle best dist to EVERY state map.
- Recipe: (1) all particles refine3D_states K=8 random, local, 15->8 A, auto its -> state_merge tree -> families (swing branch
  = states 1+3, 120.5k, 75 % pure, no truth used); (2) family K=8, maxits 8 (basis/symmetry breaker; random-label loop only
  0.002 -> 0.115 in 4 passes); (3) lrec_loop: score all vs maps (refine3D maxits=1, geom, box_crop 88, lp 8, cross-half), PCA
  of per-particle-centred scores, rebin + reconstruct + rescore 2-3x; impure family: lrec_loop2 NPC=2 (12 k-means clusters),
  motion PC = the continuous one (composition PC is skewed, 2-GMM gain 20x); (4) lrec_leaves: 10 equal-width bins of the
  motion PC (1-99 %); (5) state_merge complete linkage on leaves = chain of neighbours (all pairs z >= 10).
- Numbers: oracle 108k 0.703 -> 0.735, leaves 5.1 (4k) .. 25.3 (4k), eta2 0.54; 288k (swing30big) 0.745, 4.8 (10k) .. 25.1
  (10.7k); REALISTIC family PC2 0.50 -> 0.70, leaves 5.0 (4k, 60 % swing) .. 24.5 (96 %), eta2 0.49; ceiling true-bin maps
  0.753 (3 deg poses) / 0.772 (true poses). Per family 27k/54k/108k/288k: 0.64/0.715/0.735/0.745 -> >= ~50k per family.
- Negatives: boards even converged (15 rounds, annealed sigma; 16x1 and 7x7) corr <= 0.02, 0/16 cells reproducible, started
  at the answer they decay 0.735 -> 0.23; recursive 4x4 eta2 0.15; K=16 coordinate 0.693 = K=8; cross-half = self-half.
- Memory (user: "that's a lot of memory, is it necessary to do that many states"): MaxRSS K=8 37 GB, K=16 69, K=20 85; full-box
  scoring 327 GB OOM -> box_crop 88; request 64G for K=8 runs. 64-core jobs only fit cn001-007 -> use -c 48
  ([[batch2-csbdevel-node-classes]]).
