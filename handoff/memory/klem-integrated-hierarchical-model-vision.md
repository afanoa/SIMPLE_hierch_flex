---
name: klem-integrated-hierarchical-model-vision
description: "user 2026-10-21: the tree must be the latent model's prior, not a post hoc grouping; tree-weighted reconstruction (borrow less with tree distance), flex -> states -> focused refine -> flex per branch"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-21T05:44:07.704Z
---

2026-10-21, after I explained that the generic pipeline (f692gen: 32-state refine3D_states -> state_merge tree -> merge-up ->
coordinate sampling) is hierarchical classification without the latent model: "The whole point was to integrate hierarchical
classification with a latent variable model ... flex --> cluster states --> refine --> flex in each branch", a recursive alternative
to 32-state refinement, with "a focus mask so we have focused refinement at each node". Then: "for calculating weight distributions
for weighted reconstruction you can use the tree. If the tree organization is correct the weights will make more sense. You can
borrow less and less weights the further you go. What else could a fully integrated model give you?"

My framing (agreed direction, not yet validated): every node's map = parent + deviation; per-shell deviation powers per edge fitted to
the leaves' cross-half covariances (unbiased); leaf map = posterior mean C (C+N)^-1 x -> borrowing falls with tree distance and per
shell; the focus mask = where children differ between halves; merges = deviation consistent with zero (one criterion, cf. Bayesian
hierarchical clustering); nested soft assignment; local flex basis per node (hierarchical mixture of latent models, Bishop & Tipping
1998). Catches: junk classes need an unshrunk component; a tree repeats splits for independent changes.
First test: tree_wiener.py (routes/toy55, job tw_panel.sbatch, outputs rnasep_klem/treewiener/<set>) on the panel runs vs truth.

**Why:** the user's goal is a tree that represents the heterogeneity so merges are legitimate AND improves the maps/poses of small
states; a classification pipeline with a tree drawn afterwards misses the point.
**How to apply:** new hierarchy work should feed the tree back into reconstruction / assignment / masks; prototype each piece in
Python on the panel with truth first ([[prototype-in-python-before-fortran]]); related [[klem-bottom-up-organization-2026-10-20]].

**Results 2026-10-21 (panel runs, truth known; tree_wiener2.py, node_masks2.py, routes/toy55; outputs rnasep_klem/treewiener2, nodemasks2):**
- Borrowing per shell (posterior mean C(C+N)^-1 x) lifts every leaf vs truth: swing30hs own 7.0/10.1 A (all / 5 smallest) -> tree 4.0/3.9,
  tree + focus mask (two regions) 3.5/3.5; null_50k 9.0 -> 4.0 (no harm). Leaf-specific detail (deviation corr 40-8 A, 5 smallest):
  own 0.30, star 0.33, shuffled tree 0.36, tree 0.47, tree2 0.49, unstructured covariance emp2 0.55. prodihs (independent changes):
  star 0.05 and shuffled 0.09 ERASE small leaves (own 0.15), tree 0.40, tree2 0.43, emp2 0.44. So the organisation decides WHAT is
  borrowed; at 20 leaves the unstructured cross-half covariance is as good or slightly better than the tree-structured one.
- Truth-free organisation test: off-diagonal R^2 of the tree-structured covariance 0.85-0.90 vs shuffled 0.02-0.03.
- Focus masks per node: significance alone (negative-tail of the cross-half variance) covers 34-93 % even after per-shell regression on
  the node mean (global, reproducible differences); significant voxels holding 80 % of the variance ("sig80", focus_mask.py) = 5-25 %
  of the sphere, recall 0.8-1.0 of the true variance, 1-3x the true support (prodihs 2-10x). rel10 too small (recall 0.4-0.8).
- Hard vs M-est (user question): hard focused refinement per node, one top-2 M-estimated reconstruction of the leaves at the end;
  full soft-EM refine3D on RNase P was ~3 h/iteration (742k x 3, 80 cores) and ended at median max posterior 0.95-0.97.
- First integrated run: rflex.sbatch (flex_pca K=8 per node -> sig80 mask -> focused refine3D pcg_focus -> tree -> recurse >= 50k),
  rflex_s30a on the realistic swing30 (SNR 0.0012); baselines from lrec30/RECIPE.txt: hard K=8 oracle family eta2 0.25 (11.2-22.1 A),
  score loop on the tree-selected family eta2 0.49 (5.0-24.5 A).
- Real f692r (tree_wiener3 / node_masks2 no truth, rnasep_klem/treewiener_real, nodemasks_real): tree covariance R^2 fam1 0.41 / fam2
  0.62 (shuffled 0.08 / 0.01) -> real heterogeneity less tree-like than the panel. Borrowing keeps the arm: closed-arm occupancy
  across states corr 0.985 / 0.988 with raw, spread unchanged; closed states 16/15 0.87/0.83 -> 0.88/0.84. Member-variance masks do
  NOT localize on real data (24-100 % of the molecule) -> testing split-direction masks (node_masks3: children's mean difference).
- Split-direction masks (node_masks3, children's mean difference, sig80; rnasep_klem/nodemasks3): synthetic = the true split support
  (swing30hs 3.1 vs true 2.8 %, 1.8 vs 1.6, 7.2 vs 7.3; recall 0.85-0.95; member masks 2-5x larger). Real f692r: smaller than member
  masks (fam1 root 47 vs 70 % of the molecule, best 15 %), cover the closed arm 0.85-1.0, but stay broad (fam2 40-99 %): RNase P states
  differ in several places at once. Use member masks for a node's K-state refinement, split masks for edges / borrowing.
- rflex_s30a ROOT (400k, realistic swing30): flex_pca K=8 (14 min) = composition only (swing eta2 0.008; swing shares 61/64/46 %);
  focused refine3D (member mask 10.8 %, lp 10, 8 its, ~8 min/it) -> eta2 0.22, swing states 1 (70k, 93 % swing, 17.9 A) + 2 (49k, 71 %,
  10.1 A) = 119k at 84 % purity (recipe's random-seed root: 75 %). The 2-branch cut peeled a 20k junk state and recursed on 379k
  (wasteful) -> rflex3.sbatch DESCENT rule: peel branches < MINP as leaves, stop where the tree splits into two >= MINP families.
  Traps: this build's flex_pca rejects sigma_est/pgrp and writes the flex states as recvol_stateNN (+ _even/_odd/_unfil) in its dir.
- rflex_s30a depth 1 (rflex3 descent): r -> leaves {3 junk 20k}, r1 {4,5,6} 170k (1.5 % swing: composition), r2 {1,2,7,8} 210k.
  r2 flex states already order the swing (eta2 0.11; 79-91 % swing states at 12-18 A); focused refine -> eta2 0.375, states 8.8 A (26k,
  79 %), 13.5 A (34k, 89 %), 21.0 A (38.5k, 96 %) = beats the recipe's hard K=8 on the ORACLE family (0.25, 11.2-22.1) on a 50 %-swing
  family; score loop still ahead (0.49, 5.0-24.5). r2 -> r21 {3,5} 61k, r22 {1,2,4,7,8} 129k (swing family, depth 2 = leaves).
- Real run launched 2026-10-21 05:12: rflex_f692rf (rflex3, master rnasep.simple, ROOTROWS f692gen/seed_states.txt = 643k, mskdiam
  232), job 54035601; leaves' closed-arm readout at the end (leaf_arm.txt).
- rflex_s30a DONE 06:26 (4.5 h, 64 cores): 34 leaves; r22 (129k swing family, depth 2): flex states eta2 0.31 (flex orders the swing
  once the family is pure enough), focused refine 0.41 (8.8-22.4 A). All leaves: NMI(conformer) 0.42, purity 0.41, swing eta2 0.38,
  extent 8.8-22.4 A. Baselines: hard K=8 oracle family 0.25 (11.2-22.1), truth-seeded hard 0.47, score loop 0.49 (5.0-24.5).
  Non-swing nodes r11/r12/r21 split into ~equal flex states (no structure) -> leaves need the merge test (rflex_merge.sbatch).
  Next: coordinate stage on the continuous node (rflex_coord.sbatch: score loop from the node's focused maps, motion PC chosen truth-free
  as the smaller 2-Gaussian gain, 10 equal-width bins).
- f692gen (generic, 32-state discovery on 643k real): merge-up on the real halves kept all 32 (0 merges at every level).
- BEST SO FAR (2026-10-21 07:01): coordinate stage on rflex_s30a node r22 (rflex_coord.sbatch, 32 min: score loop from r22's focused
  maps, KL 12, NPC 2, 2 its, cross-half; motion axis chosen TRUTH-FREE = smaller 2-Gaussian gain: PC2 0.001 vs PC1 0.007) -> corr(swing)
  0.725, 10 equal-width leaves eta2 0.53, ends 4.8 A (4.4k, 60 %) .. 24.6 A (4.4k, 98 %). Recipe (tree-selected family): 0.699 / 0.49 /
  5.0-24.5; ceiling 0.753. So flex -> focused refine -> tree descent -> coordinate leaves on the continuous node beats it, all truth-free.
- f692gen DONE 08:09 (generic: board seeds -> 32-state refine3D_states on 643k -> merge-up kept 32 -> coord tree, 6 sampled branches):
  0 closed / 24 open / 8 intermediate states, closed-arm max 0.60 (cells up to 0.65, node 39). The f692r route (board families ->
  per-family refinement -> merge-up) reached 0.87 / 0.83 states and a 1.03 corner cell -> flat 32-state discovery DILUTES the closed
  state (~2 %). The recursive run rflex_f692rf is the test of whether recursion recovers it.
- rflex_s30a leaves + borrowing (tree_wiener5, truth cropped to box 88): tree R^2 0.81 (shuffled 0.05); 5 smallest leaves dev corr own
  0.60 -> tree 0.71 / tree2 0.72 (star 0.65, shuffled 0.66); merge test on the 34 leaves: 34 -> 31 (splits of non-swing families are
  conformer composition, NMI 0.42).
- OPEN: no truth-free test yet separates a continuous node from a compositional one: 2-Gaussian gains r22 0.007/0.001 vs r11
  0.003/0.000, r12 0.018/0.003; score-PCA shares and 1-D MDS fits of the node's state d2 also overlap. Coordinate bins on
  compositional nodes destroy conformer information (r11 NMI 0.01) -> deliver both per terminal node (hard states = identities,
  coordinate bins = even sampling of the leading axis) until a discriminator is found.
- HYBRID partition on s30a (hard leaves on compositional nodes + 10 coordinate leaves on the swing node r22): NMI(conformer) 0.364
  unchanged, swing eta2 0.376 -> 0.491, extent 8.8-22.4 -> 4.8-24.6 A in 10 swing states = best organisation so far.
  Path-shape test (path_shape.py: eigen shares / straightness of coordinate leaves' maps) does NOT discriminate either: r22 motion
  0.75/0.60 vs r11 composition 0.70/0.54 (the deposited conformers of the mixtures differ by motions too). Candidate left untested:
  signedness of the end-to-end difference (transport vs occupancy).
- User 2026-10-21 11:20: the 50k recursion floor is arbitrary -> OVER-SPLIT then MERGE along the tree (tree organisation is sensible).
  Plan: floor ~15k (closed-state size), depth 4, fast settings (rflex5: box 64 / nspace 2500 / 4 its / flex on <=120k / parallel
  siblings via slurm workers), then merge-up of the whole tree with the corrected test (split-decision session task_dad055d3).
  Hand-off sessions started 2026-10-21: task_dad055d3 split decisions, task_675c6104 latent integration, task_3c800824 cls_expansion
  sub-class embedding. This session keeps real RNase P (rflex_f692rf, closed-state recovery, lopsided splits, cost).
- CLOSED STATE on real f692 (2026-10-21 midday): rflex_f692rf (rflex4) 24 leaves max closed-arm 0.64 (= f692gen 0.60). Tracer (tree_tracer.py,
  f692r fam1 states 16+15 = 11.6k): kept together through the families (root 98 % -> r2, r2 86 % -> r22) but every node's 8 hard states
  spread it over ~4 states (effective number 4.0-4.5) -> fragmented at the leaves. User: tree organisation is crucial, closed states must not
  split off early -> split only where the halves support it (splitbench13 in the recursion), rare states separated late by focused steps.
- Cheap route on real nodes (cheapkids_rf*.py, view-centred flex latent PC1): r22 PC1 = arm axis without truth (closed-arm 0.11 -> 0.53
  monotone, 10 x 38k); equal-width end bins (4-6k) saturate at 0.53-0.56; r2 needs PC1 x PC2 (best cell 0.59); r21 no arm.
- Separability (closed_latent.py, LDA AUC): f692r states 16+15 in r22 whole-molecule latent 0.80 (PC1 0.76, top-2 % 3.6x); purest closed
  cell (f692r node 55 cell 25, 1.03) 0.84 / arm-focused latent (flex_pca pcg_mskfile = real-space variation mask) 0.85, PC1 0.80, 7.3x.
  Per particle the arm difference is weak -> rare states need reference aggregation (focused classification), not latent thresholds.
  Arm-focused flex PC1 children 0.08 -> 0.66 (whole-molecule 0.53). Real-space variation mask (arm_mask.py, cross-half variance of the
  arm-axis children): 16.5 % of the molecule, closed-arm region 0.74 (member masks were 24-100 %).
- splitbench13 merge on f692rf's 24 real leaves (rflex_merge_real13.sbatch, bu_merge13.py): 24 -> 20 -> 18 (r1's small leaves merge).
- Runs: focused_arm_f692rf_r22 (12-state arm-focused refine3D, seeds = equal-width arm-axis bins); rflex_f692st = rflex8.sbatch (split test
  decides every node on the real halves, decide_real.py; coarse-to-fine; 15k floor, depth 4; norm partition). f692os cancelled (fixed splits).
- CLOSED STATE RECOVERED (2026-10-21 12:45), truth-free, at the old klem's level: r22 of rflex_f692rf -> real-space variation mask (arm_mask.py on the
  arm-axis children) -> flex_pca with pcg_mskfile = that mask (focused_flex.sbatch) -> view-centred PC1 -> EQUAL-WIDTH bins: 12 bins end
  6.1k at closed-arm 0.81 (even/odd half maps 0.81/0.81, whole trajectory reproducible to 0.01-0.02); 16 bins over 0.25-99.75 % end 2.5k at
  0.89 / open 0.51 (old klem: states 0.87/0.83, cell 1.03). Whole-molecule latent ends at 0.53-0.56. Only 12 % of the end bin are the old
  closed particles (two routes, near-closed maps, different particles: per-particle membership of a rare state is uncertain at this SNR).
  Correction: the old klem region masks were DATA-DERIVED (derive_regions.py / derive_regions2.py: cross-half real-space variability of the
  board cells, components / motion clusters), not hand-drawn. Next: rflex9.sbatch (rflex_f692st, job 54041047) = focused flex at EVERY node +
  split test (user 2026-10-21).
- Tree-as-prior tests (P1 rho / P2 tree2 E-step / P3 mask units artifact / tw_units final maps): [[tree-as-prior-tests-2026-10-21]].
