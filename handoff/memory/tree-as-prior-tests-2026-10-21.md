---
name: tree-as-prior-tests-2026-10-21
description: "tree-as-prior session 2026-10-21 (s30a truth + f692u7 real): borrowed stage-B refs keep the closed state 0.87 but splitbench16 then merged it (rflex12 NOT adopted; under rflex13 focus split test both keep it); signedness and held-out likelihood do NOT tell continuous from discrete (hard wins everywhere), 1-D capture ratio does (0.55 vs 0.21-0.36); tree2 E-step references +0.014 NMI / +0.025 eta2 over raw; focus-mask 'real broader than synthetic' was a units artifact; tw_units.py = tree-weighted final unit maps"
metadata:
  node_type: memory
  type: project
  originSessionId: 34aaefbc-474e-4364-a5ac-6810fdb0bd7c
  modified: 2026-10-21T18:20:24.866Z
---

Task (user 2026-10-21): tree = the latent model's prior; three problems on the s30a truth benchmark (rflex_s30a), scripts routes/toy55/p1_*, p2_*, p3_*.
- P1 continuous vs discrete leaves (rnasep_klem/p1sign, p1hold/<node>): end-to-end difference signedness separates transport from occupancy, not
  continuous from discrete (truth r22 +0.01, r12 -0.54, r11/r21 near 0); on data biased (r22 +0.36: end bins differ in global scale/blur), a
  per-shell regression removes the bias AND r12's occupancy signal -> negative. Held-out cross-half likelihood (hard 8 / 8 equal-width bins / random
  8, one refine3D pass per half): HARD fits better on all four nodes incl. the swing node r22 (z 55-204; even on r22's swing particles) - a
  likelihood-optimised partition always wins; not a discriminator. The 1-D capture ratio rho = gain(bins)/gain(hard) over random IS truth-free and
  separates: r22 0.550 +- 0.003 vs r21 0.357, r11 0.229, r12 0.215 (one benchmark, one continuous node; a 2-conformer mix would also score high).
  Principled answer: deliver both per terminal node - hard states = the likelihood model (assignment, borrowing), coordinate bins = sampling of the
  leading axis (ends 4.8-24.6 vs 8.8-22.4 A), with rho reported. Self-inclusion ([[refine3d-lp-set-scores-merged-map-trap]]) measured negligible.
- P2 tree weights inside refinement (rnasep_klem/p2estep/<arm>, p2_arm.sbatch: refine3D maxits=1 x2 over 400k x 34 leaves, references = tree_wiener5
  estimates): after 2 its NMI raw 0.424 / own 0.427 / tree2 0.438 (start 0.416), swing eta2 0.367 / 0.371 / 0.392, pose median 2.71 / 2.68 / 2.62 deg
  (start 2.63: no pose room), new leaves' dev corr all/5 smallest raw-map 0.531/0.392 vs tree2 arm 0.569/0.432, after final borrowing 0.680/0.575 vs
  0.751/0.707. Gain = borrowing, not shrinkage (own = raw). Cost: borrowing ~30 s/it vs 44-45 min refine3D.
- P3 focus masks (rnasep_klem/p3masks, p3blur, p3junk, p3views, p3poserr): node_masks3 reported SYNTHETIC masks per sphere (3.1 %) and REAL per
  molecule (15-95 %); the molecule is ~3 % of the sphere -> in matched units real masks are NOT broader (member masks: panels 50-100 % of the
  molecule, s30a realistic flex nodes 97-100 % with true support 18-70 %, real f692rf 87-100 %, f692r fam1 32-86 %). At realistic SNR the member
  masks over-cover in BOTH (s30a r22: 99 % vs true 17.6 %, 74 % of the reproducible variance off the true support). Not the cause: view
  distributions, junk, local blur/scale, rigid misregistration (0.8-1.0 deg fitted, no change), per-component masks. The excess is low resolution:
  s30a r22 band masks 40-20 A 83 %, 20-12 A 26 % (~true); real f692rf r 20-12 A 27 % but loses the closed arm (0.17). So the rflex focused
  refinements were effectively unfocused; focus only where a mask is local in MOLECULE units.
- Final maps (relayed task, rnasep_klem/tw_units/f692rf_merge13): routes/toy55/tw_units.py + tw_units.sbatch (env MERGE, RT): unit half-maps
  (count-weighted leaves) -> unit tree (state_merge --list) -> raw/own/tree2/emp2 at box 160 lp 10 unsharpened -> view_rflex/<RT>_merged_tree2,
  _emp2 (halves/), choice by tree R^2 >= 0.7. merge13: R^2 0.741 -> tree2; closed-arm corr with current 0.990 (emp2 1.000), spread 0.166 vs 0.180,
  calls unchanged 0/5/9; units >= 20k so 10 A barely changes (own weight 0.91 at 20 A, 0.55 at 10 A).
- (b) borrowed references INSIDE a node's focused refinement (routes/toy55/tw_stageB.sbatch / tw_stageB2.sbatch: RT ND ARMS=raw|auto
  START=flex|stageB; refine3D maxits=1 x4 with the rflex4 node keys incl. pcg_focus, refs rebuilt each iteration; rnasep_klem/tw_stageB/s30a_<node>):
  the one-iteration raw loop reproduces the native single call (r22 eta2 0.403 vs 0.41). Node trees of 8 siblings fit poorly (R^2 0.24-0.50) ->
  auto = emp2 always. After 4 its raw -> emp2: r22 NMI 0.128 -> 0.139, swing eta2 0.403 -> 0.432, leaf dev raw-map 0.482/0.432 -> 0.558/0.506
  (after borrowing 0.670/0.635 -> 0.744/0.715); r21 NMI 0.155 -> 0.172, purity 0.484 -> 0.493, dev 0.321/0.217 -> 0.365/0.270; r12 NMI 0.322 ->
  0.337, purity 0.523 -> 0.534, dev 0.603/0.565 -> 0.631/0.603. Iterations 3-5 min; borrowing + restart ~1 min. Proposal: rflex11 = rflex10 with
  BORROW=auto stage B as single-iteration calls (do not edit rflex10); Fortran port = posterior-mean references in the master's volume assembly.
- Real-data minority safety (rflex_f692u7 node r, tw_stageB3.sbatch + node_readout2.sbatch, rnasep_klem/tw_stageB/f692u7_r): borrowed stage B
  (emp2; node tree R^2 0.15) keeps the closed state: state 6 8468 ptcls closed-arm 0.87, old closed 15+16 19.8 % of it (raw loop 8745 / 0.87 /
  19.7 %, native 8821 / 0.87 / 20.0 %); MSKDIAM 232 + local parts reproduce the native run.
- routes/toy55/rflex12.sbatch = rflex11 (another session's: rflex10 + splitbench16 + DECIDE) + BORROW auto|emp2|tree2|none (stage B as ITB single-
  iteration calls, refs rebuilt each iteration, fstageB/borrow_it<k>, it_done/stageB.done resume). END-TO-END node test (rflex12_seedtest.sbatch:
  new runs rflex_f692u7b / f692u7n seeded from f692u7's root up to stage A, MAXD 0): stage B as above, BUT splitbench16 then merged the closed
  state 6 with state 4 in the BORROW run (leaf 4+6 15.9k at closed-arm 0.68) while BORROW=none kept it alone (8.8k at 0.87). The beyond-EM re-test
  runs at its noise floor there (pair statistics 0.00-0.05 vs null max 0.05), and borrowed states are less distinguishable over the full band
  (pair statistics 0.63-0.79 vs 0.68-0.84) -> NOT adopted for node refinement yet; borrowing stays for final maps (tw_units) until the decision step
  is robust to it.
- 2026-10-21 16:58 (other session): f692st resumed on routes/toy55/rflex13.sbatch = rflex12 + split test restricted to the node focus mask (--focus rsmask.mrc --rb_margin 4: keeps all 8 f692u7 arm states, off the noise floor) + falloff-prior node masks (arm_mask2.py --frac 0.2 --dil 3 --falloff 6), run with BORROW=none (job 54048912); borrowing stays for final maps (tw_units.sbatch on the finished merge). Untested: BORROW=auto under the focus-restricted split test.
- 2026-10-21 17:25 rflex13 seeded node test (routes/toy55/rflex13_seedtest.sbatch, runs rflex_f692u7b13 / f692u7n13 from f692u7 stage A, old hard rsmask kept): with the FOCUS-MASK split test both BORROW=none and auto keep all 8 states (K* 8) -> closed leaf state 6 none 8816 / 0.87 / 20.0 % old-closed (37.6 % of the node closed), auto 8473 / 0.87 / 19.7 % (35.7 %); all leaves within 0.01. So the rflex12 merge was the unfocused split test at its noise floor; under rflex13 borrowing is safe for the closed state at this node (3.9 % smaller closed leaf, same purity), no real-data gain shown. Untested: falloff masks + borrowing.
