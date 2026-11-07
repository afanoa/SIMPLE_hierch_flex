---
name: klem-bottom-up-organization-2026-10-20
description: "fast-loop boards are flat (maxr ~1/49, 95 % argmax change) and miss within-family motion; K-state refine3D_states inside a family resolves the swing (random seeds as good); bottom-up tree from cross-half difference energy of leaf half-maps (state_merge.py) merges only neighbours; vol1..vol20 caps K at 20"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-20T23:52:00.256Z
---

Evening 2026-10-20, toward the goal in [[klem-organization-and-merge-up-2026-10-20]] (tree where any merge is legitimate).
- Fast-loop klem boards (5 rounds, 100k synthetic, SNR 0.005) are flat per particle: median maxr 0.02-0.05 (1/49 = 0.020),
  argmax changes 95 % every round; in the swing family (ff3_swing30hs fam6) region-board cells explain 0-9 % of the true
  swing variance (only region_5 any). Factor bins 10.1-18.5 A of 0-30 were this weak signal, not a sharpening problem.
- K-state refine3D_states INSIDE the family (fam_classify.sbatch, pose_policy local, 15->8 A, 8 its, ~12 min on 24k):
  8 states ordered along the swing, means 4.6..22.7 A (bins seeds) / 6.0..21.4 A (RANDOM seeds), within-state sd 3-5 A
  vs 6.9 family -> the boards are not needed for this step. The family's own swing is peaked (sim ENDW 0.35) and its
  extremes sit in other K=8 composition states (low end with apo = the closed conformer, high end 24.7 A with apo).
- Bottom-up organisation: state_merge.py (routes/toy55; state_merge.sbatch) on leaf even/odd maps: num = Re<M(E_a-E_b),
  M(O_a-O_b)> over 40-10 A (unbiased difference energy), z = num/sd (null verified ~N(0,1) on toy maps), d2 = num /
  sqrt(signal_a signal_b) (size-independent). z is useless as a merge rule (every pair of 5k+-particle states z >= 50);
  order merges by d2. On the 8 swing states every merge joins swing neighbours, ends kept to the 4-state cut
  (4.6 | 11.8 | 17.3 | 22.7 A). Centroid linkage CHAINS on real data (f692r 218 leaves: one 69k cluster at the 84 cut);
  average / complete (--link) keep clusters balanced.
- f692r old tree measured by closed-arm eta2 (share of leaf closed-arm variance explained): zoom nodes (63) 0.90,
  region boards (11) 0.12, families (2) 0.11, merge-up 84 0.92 -> the upper levels do not organise the arm.
- f692r bottom-up on the 218 zoom leaves (consolidated_12_fullw_zoom halves, 1.5 min, 8 cores; f692r/smerge_zoom218_<link>,
  cutN/{groups,states}.txt, tree.nwk): closed-arm eta2 complete 84/48/24/12/6 = 0.92/0.83/0.78/0.57/0.55, average
  0.91/0.85/0.70/0.67/0.24 vs the old tree's region boards (11) 0.12. Complete linkage keeps the closed leaves 22 (0.931) and
  4 (0.921) in separate states to the 12 cut and joins them in one closed branch (30 leaves, 74.5k, mean 0.59) at 6.
  Figure: routes/toy55/plot_tree.py (merge_tree.txt + arm_states.txt). Reconstruction of the complete cut48: mu_rec.sbatch RTAG=f692r
  CDIR=smerge_zoom218_complete/cut48 CUT=bu48.
- Synthetic K=20 refine3D_states (fast_comp KSEED=20, ~36 min on 100k) -> state_merge: swing30hs purity 0.96, tree keeps
  class purity 0.93-0.96 to 8 states, the 4 swing states form one branch merged in swing order; apo-holo pairs merge early
  (holo = apo conformer + bound part: the tree orders by SIZE of difference, not type -> composition is not automatically on
  top). prodihs (swing x breathing grid, hs): grid states merge as (k, j) neighbours, purity 0.99 to 6 states.
- REALISTIC SNR (0.0012): within-family K=8 on 30k swing particles resolves nothing (means 13.8-17.6 A); oracle family of all
  108k swing particles K=4: 11.2..20.5 A, K=8: 11.2..22.1 A (only the far-swung end separates, sd 5.7-7.2) -> at realistic SNR the leaves (classification) are the limit, not the tree.
- Map-level z is INFLATED for board-derived leaves (labels used both halves): f692r leaf pairs min z 8.9, nothing merges ->
  merge DECISIONS need the merge-up's label-shuffle null (mu_family: all pairs within a family); d2 only orders the tree.
- Deliverable: merge-up 84 (f692r/mergeup/consolidated_family, tested merges) + bottom-up tree over them
  (f692r/smerge_mu84_complete): closed states 15 (0.861) and 16 (0.896) are nearest neighbours (merge 30), nested in the
  intermediate-arm branch; closed-arm eta2 48/24/12 = 0.95/0.86/0.79. Budget cut (--budget) folds small intermediate leaves
  INTO closed leaves (4 -> 0.76) -> extremes need particle re-assignment, not map-level folding.
- GENERAL PIPELINE (user chose option a = all-particle discovery + cross-half decisions, 2026-10-20 evening): panel_run.sbatch
  (K-state refine3D_states, random seeds, all particles) -> state_merge.py tree -> bu_merge.sbatch (tree cut N/2..1, per group
  splitbench9 on even/odd rows of ONE stack (cand files per half, no stack copies; eo from the refinement's OUTPUT project - input
  copies have eo=0), --pair identity --xem 3 --strata 1) -> scorecard.py. NULL CONTROL null_50k (1 conformer + 2 junk, K=8):
  free or view x defocus stratified shuffle null keeps 7 states (conformer splits reproduce across halves: shared noise from
  all-particle discovery; pose err/corr/shift/contrast/scale identical); with the cross-half EM (xem 3) -> 3 states = truth
  (conformer 99.9 % pure). The cross-half EM is mandatory in the decisions. two_50k (2 conformers + 2 junk, K=8) -> 4 = truth.
  swing30hs board-seeded K=20 leaves -> 14 states (class purity 0.96 kept, NMI 0.65 -> 0.69, swing 3 states 6.1/14.9/21.7 =
  ends kept); RANDOM-seed K=20 discovery is worse (purity 0.93, swing 2 states) -> panel_run2.sbatch seeds from the whole
  board's cells (spectral K clusters; board reused from ff2_<TAG>/whole/board or built). Discovery under-segments
  conformers (16 apo/holo -> ~9 states): needs 2-level discovery or K > 20 (handoff).
- SAMPLING (general): state_axes_tree.py = tree over the FINAL states (half-maps = particle-weighted leaf averages) +
  cross-half principal axes at every node with >= 3 states (leading axes to 90 %, <= 3), states ordered along each axis.
  swing30hs: the swing branch's axis 1 (79 %) orders 6.1 -> 14.9 -> 21.7 A (rho +1.00), ends = the extremes.
- COORDINATE METHOD folded in as stage 4 (from the realistic-SNR session's lrec_loop/lrec_pc: cross-half scoring of particles
  against K state maps (refine3D 1 it, prob tables), PC1 of per-particle centred scores = coordinate, rebin, iterate):
  copies coord_pc.py (generic --codes truth + even-width bins 2nd..98th pct, small bins folded inward), coord_loop.sbatch
  (OUT/TAG params, XHALF 1), coord_nodes.sbatch (per node: >= 3 states start from the states, a single state from a random
  8-split; auto = final states >= MINP), coord_tree.sbatch + coord_tree.py (ONE cross-half scoring pass of all particles vs
  all final-state maps, nspace 5000 box 88 -> PC1 per tree node from the score sub-block, no re-scoring).
  swing30hs swing branch (3 merged states 6.1/14.9/21.7 A, 28.9k): corr 0.88, 8 even bins 3.6..25.8 A, 2 % tails 2.9/26.8,
  eta2 0.62 -> 0.78, 10 min. The merge-up's cross-half EM drains continua (prodchs 13 grid leaves -> 3 states); the coordinate
  step restores the extent -> design: merge-up (discrete, conservative) then coordinate (continuous within).
- SAMPLING RULE (coord_tree.py; v2 = routes/toy55/coord_tree2.py): per node of the final-state tree the fewest axes (<= 3) carrying
  70 % of the score variance; samplable if no axis has a deep valley (24-bin histogram valley < 0.6 of the smaller flanking peak);
  even grid 8 / 4x4 / 3x3x3, small cells folded into the nearest kept cell; sample the maximal samplable nodes. 1-axis-only rule
  (share >= 0.65) picked only state pairs on real data -> real branches are multi-axis. Synthetic: swing branch 4.0..25.5 A
  (one scoring pass ~= iterated loop 3.6..25.8), coupled grid branch 1.0..11.1 of 0..15; within-state random-split coordinate
  weak (breathing corr 0.31); realistic 100k swing fails (corr 0.04, 26.7k swing particles < ~50k support).
- f692r REAL (coord_tree_real.sbatch per family, focus3 snapshot for > vol20, mskdiam 232; recon_nodes.sbatch for cells; maps in
  f692r/coordtree_fam{1,2}/nodes[2]/nodeNNN/maps/state_NN.mrc, halves in halves/): fam1 (30 states, 240k) 1-axis branches = arm
  segments (e.g. node 32: closed-arm 0.56 -> 0.20, open 0.92 -> 1.42); closed pair 15+16 coordinate 0.83..0.94. Multi-axis fam1
  -> nodes 54/55/56 (3 axes); node 55 (closed pair + near-closed + intermediates, 12 states, 106k, 27 cells) spans closed-arm
  0.18..1.03; cell 25 (corner +1.4/+1.3/-1.4 sd, 2.5k) = closed-arm 1.03 / open 0.55, the most closed RNase P map so far.
  node 54: closed-site and open-site occupancy vary independently (4 corners), not a rigid two-site swing. fam2 (54 states)
  multi-axis -> 9 branches, arm 0.07..0.55 (no closed states in fam2).
- Focused PCG classification (refine3D nstates 8, rec_backend=pcg pcg_mskfile=region) COLLAPSES: populations 548..7527 ->
  9..10617 in one iteration, then "distributed PCG requires both halfsets"; the unmasked refine3D control stays stable.
- SIMPLE: only vol1..vol20 exist as command-line keys -> refine3D_states K <= 20 (K=24 dies in prob_tab_neigh:
  "vol21 argument is not allowed"). pcg_mskfile is a HARD support (out-of-mask density removed; the fixed-background
  focused mode x = x_fix + P delta is deferred, not built), so a region mask gives region-only references.
