---
name: node-type-test-disjoint-band-2026-10-21
description: "truth-free per-node test (relayed task): STOP works = z of the cross-half reproducibility of a node's groups in a band the grouping never used (cheap PC1 children 15-11 A, hard states 10-8 A; 8 Fourier-sector SE; go if max z >= 4); continuous vs compositional NOT decidable from map geometry even on truth maps; hard states capture more reproducible heterogeneity everywhere -> default hard states"
metadata:
  node_type: memory
  type: project
  originSessionId: 6d5934fe-4a66-44cc-9cf7-17b64c67240a
  modified: 2026-10-21T19:34:13.924Z
---

2026-10-21, task relayed by the 2D-to-3D session for the user (node type -> recursion: continuous = coordinate bins, compositional =
hard states, nothing = stop). Scripts routes/toy55: nodetype_stop.py/.sbatch (children), nodetype_stop2.py/.sbatch (hard states),
nodetype_stats.py (children path geometry), nodetype_states.py (hard-state geometry), nodetype_compare.py (H = member-weighted
cross-half deviation power of two organisations), cheapkids2.sbatch (truth optional, MSKDIAM), nodetype_devsetup.sbatch (junk-free
pseudo-nodes rc of null_50k / two_50k). Outputs rnasep_klem/nodetype_*.out, cheapkids_{s30a,root_null_50k,root_two_50k,f692rf}_<node>.

- STOP (works): assignment band never overlaps the measurement band (flex latent <= 15 A -> measure 15-11 A; refine lp 10 -> 10-8 A);
  per-shell regression of every group on the weighted mean first. Dev z: nulls 1.0-1.8 (null_50k junk-free, null root, r11);
  heterogeneous max(children, hard) 4.7-27.8; two_50k junk-free 10.7. FAILS when PC1 is a junk/nuisance axis: two_50k root with
  junk z 2.9, real f692rf root children z 1.2 while its hard states z 15.6 -> use both routes, run after junk removal. Cannot certify
  heterogeneity that lives only coarser than the free band (false stop risk). The biased (assignment) band is useless: r11 R 0.44.
- v3 (traps from the split-decision session; nodetype_stop3.py/.sbatch): bands 2 shells past the assignment edge (children 13-9.5 A,
  hard 9-7 A) and a per-group 3-D shift fit against the other half's consensus. Shift fits are 0.02-0.27 px and move z little, although
  the PC1 children DO carry a monotone pose-error trend at heterogeneous nodes (r22: rotation-error means +-0.6 deg end to end, shift
  0.3 px; r11 none; pose_err_groups.py). Two routes ~independent (different bands) -> Stouffer z_comb = (z_children + z_hard)/sqrt2,
  GO if z_comb >= 3 (single route: z >= 3). Dev: true nulls z_comb <= 1.4 (null_50k junk-free 1.4, r11 1.3, s30split2 r1 -1.6, r2 0.5);
  go s30a r12 7.8, r21 6.3, r22 3.9, two_50k junk-free 5.9, s30split2 r3 4.2, r5 4.3; FALSE STOP s30split2 r4 0.8 (holo 101/103/107,
  refined NMI 0.12). s30split2 r1 (93 % swing over 15.5-19.9 A) refined swing eta2 only 0.05 -> stop is right; split test K* = 1 too.
- SIDE BY SIDE with the split session's joint K* (rflex6split_decide over splitbench16, consensus + node poses; rflex_s30split5):
  agree on 10 of 11 shared nodes; the miss is s30split2 r4 (their K* 2: 47k + 8.9k state, EM-children NMI 0.17/0.14). r3 K* 2, r5 K* 3,
  r1/r2 K* 1. The gate's global statistic has little power for ONE small distinct state; their pairwise equal-N test catches it.
  The cheap children route alone cannot be a PRE-refinement filter: it stops the real f692rf root (z 1.4, PC1 nuisance), s30split2 r4
  (-0.5) and r5 (2.6). So the gate is a post-refinement second opinion; K* from the split test should carry the split decision.
- REAL f692rf calls (v3 z_comb): r 6.9, r1 3.1 (borderline), r2 8.9, r21 11.0, r22 5.7. Earlier (v2, max z >= 4): r go (hard z 15.6), r1 STOP (children 2.1 / hard 3.7; r1 = root state 7, flex gave equal states), r2 go
  (6.0 / 10.7), r21 go (7.2 / 7.5), r22 go (4.0 / 6.0).
- PER-STATE GATE (user ask; nodetype_perstate.py/.sbatch, out nodetype_perstate.54045687.out, synthetic perstate_synth.py): every
  hard state's halves vs the OTHER half's leave-one-out mean, a Newton-refined 3-D shift fit of every state first (the parabolic
  peak of a band-pass correlation is biased the same way in both halves: it flagged a 0.6 px shifted synthetic state), per-shell
  regression, 32 equal-solid-angle sectors, t(31), Bonferroni over K (z >= 4.03 at K 8). Synthetic: nulls max z <= 2.6, blur/shift
  nuisances clean, one 6 % state z 5-7 where the global z is 1.0-1.7. Dev 9-7 / 9-6 A: nulls r11 2.3/1.4, s30split2 r1 2.3/1.9,
  r2 1.6/1.9, null_50k 3.9 (a junk state)/2.5 = all STOP; GO r12 9.7, r 43, r1 11.5, r2 20.6; only at 9-6 r22 5.1 (swing-end state),
  s30split2 r3 4.1 (= the split session's 9.3k state 4), two_50k 4.8 (= the conformer-18 state); MISSED r21 3.0/3.5 (diffuse: all
  states mixtures), s30split2 r5 2.9/2.7 and r4 1.7/2.0: r4's top state IS the split session's 8.9k state 7 but the free band holds
  no signal (global 1.7) -> the r4 miss is the BAND (states assigned at lp 10 leave 9-6 A; the split test labels at <= 20 A and tests
  17-10 A), not the aggregation. Pure-noise junk states are not flagged (z <= 1, fitted shifts 5-13 px). Put into the v3 Stouffer
  instead of the global hard z: same calls on all 9 nodes (r5 3.2/3.0). Real f692rf: r 15.0 (8/8 flagged), r1 2.0 STOP, r2 16.9
  (8/8), r21 5.3 (states 4, 8; 3.8 at 9-6 vs global 9.5 = diffuse), r22 6.8 (6 flagged). Use it to NAME the distinct states, not
  as a replacement for z_comb.
- LP-20 STATES FIX r4 (user ask; lp20ref.sbatch -> rnasep_klem/lp20_s30split2_<node>; nodetype_perstate2.py adds the truth deviation
  share tdev and the overlap with another labelling): the node's focused refinement repeated with ONLY lp 20 (matching 19.52 A, 13-15
  min @ 64 cores) frees 17-10 A. Per-state GO at all three: r3 8.5 (states 4, 6; at 13-10 A, past the flex <= 15 A start labels too,
  5.8, state 4), r4 7.8 (7, 8, 1; 13-10 A 5.4: 7, 8), r5 8.4 (5, 1, 6, 7; 13-10 A 7.0: 5, 1); the global z also GOes (5.1/6.0/6.4, was
  2.8/1.7/3.5 at 9-7 A on the lp-10 states). Top states = the split session's: r3 state 4 (75 % from lp-10 state 4, 60 % conformer 10),
  r4 state 7 (77 % from lp-10 state 7, 67 % conformer 107). Flagged = the conformer-enriched states (51-67 % top); z vs truth deviation
  share Spearman 0.71/0.98/0.74. Truth deviation is LARGER at 9-7 A (r4 state 7 0.15 vs 0.09): the 9-7 A miss was data SNR. lp 20
  costs nothing on truth (NMI r3 0.09 -> 0.10, r4 0.12 -> 0.16, r5 0.18 -> 0.17). NOT yet null-calibrated at lp 20 (r11, s30split2 r1/r2).
- CONTINUOUS vs COMPOSITIONAL (not achieved): on noise-free TRUTH maps the children paths are equally linear (r22 0.90 vs r12 0.90,
  r21 0.91), neighbour profiles U-shaped for r22 and r21, hard-state 1-D share r22 0.66 vs r/r1/r21 0.69-0.74, H_hard/H_bins r22 1.50 vs
  1.32-2.47, H_hard/H_scoreloop r22 1.56. Groups along PC1 are mixtures either way and r22 itself is 79 % swing + 21 % conformers.
- DECISION: hard states capture more reproducible heterogeneity than bins at every heterogeneous node (dev truth 1.3-5.0x, real data
  3.4-118x) and costs are asymmetric (bins at compositional nodes NMI 0.34 -> 0.10 r12, 0.15 -> 0.08 r21; hard states at the swing node
  eta2 0.53 -> 0.41) -> hard states at every resolvable node; decide continuity per hard-state FAMILY (score loop + coordinate bins of a
  family, as rflex_coord at r22).
Related: [[subclass-node-embedding-2026-10-21]], [[klem-integrated-hierarchical-model-vision]].
