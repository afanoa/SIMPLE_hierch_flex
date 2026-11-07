---
name: cahra-13991-pose-entanglement-run-2026-10-23
description: "CAHRA 13991 (CALHM2 C11, 30k simulated, 18.5 % open): the sim's two models differ by +16.2 deg in frame, so the consensus superposes the conformations and the open poses are correct (1.3 deg median in their own frame); C1 trees fail, and C11 plus the contrast fix (m8c26) separate open at the root; work dir batch2 afan/cahra_13991"
metadata:
  node_type: memory
  type: project
  originSessionId: 5f12aee4-62b6-4518-9284-7a4df90b8e47
  modified: 2026-10-25T03:12:24.179Z
---

Run 2026-10-23 at the user's request, coordinated with the "2D to 3D heterogeneous state translation" session.
Work dir: batch2 /mnt/beegfs/elmlund/afan/cahra_13991. Data: EMPIAR_downloads/13991/data.
- Data: Particles/pose_entangle.mrcs (30k, box 256, 1.2 A/px, float32). pose_entangle.star has defocus but no angles. source.star
  is the truth: conformation 0 = closed 6UIW (24,450), 1 = open 6UIV (5,550). Nearly all views are side views.
- Import: stk + plaintexttab (DefocusU V Angle, dfunit=A angastunit=degrees), avoiding the STAR traps.
- TRAP: m8c15 has no `abinitio3D`. Ab initio there is `solve3D` (also in ~/SIMPLE_proc).
- TRAP: EBI serves ~1 MB/s per connection. Use 12 parallel `curl -r` ranges (16 MB/s).
- Consensus (solve3D c11, then refine3D_auto c11): closed poses 0.7 deg median. Open poses sit at a consistent +16.1 deg about the
  C11 axis in the global (closed) truth frame. That is NOT a pose error: the sim's open and closed models are themselves rotated
  +16.2 deg in frame (truth-frame map docking), so the consensus correctly superposes the conformations. In their own frame the open
  poses are 1.3 deg median. The deposited 6UIV/6UIW frames differ by 9.4 deg, the sim's step7 models by 16.2 deg.
- hetero_tree in C1 (m8c15, V3 gate, minpop 3000, maxd 3), seeds 1/2: every leaf at the base open rate, NMI 0. The replicates
  disagree (one leaf vs a 14k/12k split). The per-site pole (6,476, halves z 22) is not the open state (14.8 % open).
- Oracle at the consensus poses, true labels: C1 maps reach only 11.8 A (closed) / 18 A (open); C11 maps reach 2.40 / 2.95 A.
  Open-vs-closed band correlation in C11 is 0.94-0.72, vs a random split's 1.00-0.94. The difference lies 10-30 A from the axis
  at mid-height (pore N-terminal helices). In C1 the signal is only at 30-15 A.
- STAGE SCORING (stage_score.py; root node, both seeds): C1 flex DOES find open (z4 AUC 0.83, CV-LDA on z1..z10 0.88; flex label
  5 = ~3k at 74 % open). stA0 keeps part of it (best state 50 %). Then stageA/stageB (K=2) and the seeds split by the truth
  per-particle SNR (AUC 0.04) at 18 % open on both sides. Each CAHRA particle has its own SNR (_rlnsnr 0.003 to >0.02).
  So the loss is per-particle amplitude in the state assignment/seeding plus the stA0 -> stage A grouping, not detection.
- Truth-free confirmation (stage_score2.py): the junk-screen amplitude has Spearman 0.84 with the truth SNR, and it predicts the
  stageA/B and seed labels (AUC 0.14-0.17, medians 0.88 vs 1.10). The junk screen flags 3.3k GENUINE low-SNR particles (11 %; the
  set has no junk). Port diagnosis: stage A/B (euclid) and the site EM have no per-particle amplitude; seeds come only from flex
  z1/z3 tails (site_dirs '1,3'), so the open z4 direction never seeds. The fixes need the user's OK.
- m8c15 + pgrp=c11: flex/refine/reconstruct use c11, but state_split, landscape_sample and state_mergeup run c1. The site EM,
  junk screen and landscape cell_of use no symmetry in any build; symmetrising them is a code change that needs the user's OK.
- C11 trees (job 54119722; pgrp=c11, split test still c1): 2 leaves, 13.3k pure closed vs 16.7k at 33 % open (98.5 % of open),
  NMI 0.19, replicates NMI 0.70, 0 junk. Open is FOUND 4x: root flex label 4.2k at 82 %, root stA0 ~5k at 76-79 %, r1 flex 3-4k at
  76-82 %, r1 stA0 4-5.5k at 66-80 %. It is LOST at the root 6->2 grouping (which follows amplitude: AUC 0.8) and at r1's split
  decision (K=1). The site poles reach only 49-58 %.
- site_dirs diagnostic (C1 root, site_dirs=1..10 on fflex AND flex, sitedirs_54126468.out): all 20 endpoints, including flex z4hi
  (latent AUC 0.83), converge to the same 2 fixed points (6k at 19.7 % open, 15k at 17.2 %): the amplitude split. Seeding more
  directions does not help until the site EM has a per-particle amplitude. Cost ~40 s/direction on 32 threads (~21 core-min).
- Fix builds: m8c21 = pole rule; m8c23 = per-particle contrast (site EM + stage A/B) + symmetry in the site EM / junk / landscape;
  m8c24 = identity-test junk rule. On CAHRA (m8c21 C1 root replay, stree_port/junkval_cahra_c1_m24) the m8c24 rule passes: 0
  flagged, the faint tail kept (3 amplitude components, truth SNR medians 0.0043 / 0.0088 / 0.0207); the old rule flagged 4,849
  genuine. The reruns are tree_m8c21_c1_s*, tree_m8c23_c11_s*; tree_rerun3.sbatch + stage_score3.py add the contrast checks.
- TRAP: the m8c24 identity-test junk rule passes in C1 but FAILS in C11 (m8c26 root, both seeds): it flags the faint component
  (amplitude 0.70, ~9,950 particles, 30 %) as "a different structure" at z 18 although that map correlates 0.953 with the genuine one.
  C11's 11x averaging makes small systematic (blur/pose) differences significant, and the rule has no effect-size floor. Validate
  junk rules in the point group the tree will run, not C1.
  Diagnosis (junk_diag.py): flagged = the faint tail (truth SNR median 0.0045 vs 0.0133; 20.2 % open vs 17.8 %, not the conformation).
  Poses are worse because of faintness (median 10 deg vs 0.6; within SNR bins 13.2 vs 12.4, 5.7 vs 4.5). After per-shell scaling the
  residual difference is ~0.03 % of the map power, flat across 80-10 A; corr with the true open-closed difference is 0.11. The test
  needs an effect-size floor.
- m8c26 C11 ROOT (the contrast fix works): stage B = 4.7k at 74-75 % open vs 16.2k at 1.4 % (NMI 0.56-0.58; ~94 % recall of the
  kept open). Stage A/B amplitude AUC 0.41 (was 0.8), truth-SNR AUC 0.45; fitted contrast Spearman 0.89 with truth SNR. The leaves
  come from m8c27 (decide fix, old amplitude junk rule).
- Twist (pose_stage.py): open particles stay at +16.1 deg about the axis (median 15.7, p90 16.4) through stA0/A/B, even inside the
  75 %-pure open stage-B state. It is a tight, self-consistent per-state frame (about half a C11 subunit); refinement cannot move it.
  RESOLVED (dock_perstate.py, truth-frame maps in cahra_13991/gtframe/gtp): the two SIMULATED MODELS differ by +16.2 deg about the
  axis at best superposition (= the open poses' offset +16.0). The "twist" is a frame convention; the open poses are correct
  (own-frame error 1.3 deg median). The consensus superposes the conformations. The closed contaminants of the open state are fainter
  (truth SNR median 0.0083 vs 0.0138) and partly misaligned (p90 18 deg).
- m8c27 C11 (decide fix, old junk rule; job 54140378): the FIRST open leaves, 2,990 at 85 % open (s1, recall 46 %) and 2,174 at 87 %
  (s2, 34 %); closed leaves 98.7-99.2 % pure; replicates Jaccard 0.62 on the open sets. Recall is capped by the 9k root junk
  (1,825 open lost): the open leaf gets ~68 % of the open that survived the screen. m8c30 (empty-only junk, 0 flagged) is job
  54148398. Scorers: leaf_pop.py, tree_rerun4.sbatch, pop_score.sbatch.
- m8c30 C11 (empty-only junk = 0 flagged; leaf bug present): root stage B holds 91 % of open, but with 3.3k faint, misaligned
  closed (purity 60 %). Leaves: open-majority r12 3.6k at 58-61 % (recall 37-41 %; with the leaf fix ~59 % / ~56 %). TRAP at r2:
  stage B separated 2,094 at 74 % open (NMI 0.63), but state_split returned K*=1. CAUSE (port): a fixed 1000-particle floor per half
  (EM_MIN) in the cross-half EM; the minority group had 998 / 986 per half, so it was dropped. My "band too coarse" guess was WRONG.
  The floor is removed in m8c38em. m8c36 (leaf fix) is job 54172221 (tree_rerun5.sbatch + leaf_cause.py). m8c36 can crash with
  site_layer=yes (fsite_* collision, EOF in site_seeds.txt); fixed in m8c38. A hetero_tree resume must name the input project as
  copied into 1_hetero_tree (p.simple), not root.simple.
- ORACLE CEILING (oracle_ceiling.py; true maps, per-particle best-explaining map, 40-6 A, slice A^T): true poses accuracy 0.982,
  NMI 0.83 (faintest quintile still 0.95); consensus poses accuracy 0.929, open precision 0.77, recall 0.88, NMI 0.53 (Q1 0.80). The
  m8c36 tree (open leaf 70 %, recall 0.51, NMI 0.11-0.15) is far below the consensus-pose ceiling: the gap is the tree, not faint
  particles.
- m8c36 C11 FINAL (leaf fix): all 30,000 placed (0 junk, 0 unseeded). Open leaf 4.0-4.1k at 69-70 %, recall 0.51; replicates
  Jaccard 0.97. Root stage B already holds 90 % of open (68 % pure, at the consensus-ceiling recall); the descent loses ~40 %:
  root split groups send 1.5k open-state particles into r1 (s1), and stA0 drops a small 1.2k open fflex label (s2 r1).
  Next is m8c39 (EM floor removed; job 54180051, tree_rerun6.sbatch).
- m8c39 C11 (EM floor removed; job 54180051) = BEST: all 30,000 placed; open-majority pool precision 0.73-0.75, recall 0.71-0.75,
  binary accuracy 0.90, NMI 0.39 (consensus-pose ceiling 0.93 / 0.77 / 0.88 / 0.53). The floor fix recovered r1_s1 1,718 at 78.5 %
  (s1) and r12_s1 1,176 at 94 % (s2). The remaining loss: r2's split routes ~950 open-state members into r21, where nothing
  separates; that is the port's SIMPLE_SPLIT_ROUTE=best fix. Scorers: split_trace.py, early39.sbatch, tree_rerun7.sbatch.
  Cluster host is now batch.ncifcrf.gov (2026-10-24). Arm scoring: score_arm.sbatch (TAG, T1, T2 = .../1_hetero_tree dirs) ->
  arm_<TAG>.txt; compare_arms.py GT TAG:T1:T2 ... gives one table vs both ceilings (m8c39 baseline in arm_m8c39.txt).
- ARMS M / MF (m8c52t = m8c50 + SIMPLE_SPLIT_CONTRAST=maps; MF adds flex_min_neff=1000): the contrast-aware root split sends 4,758
  of 5,550 open into r2 (6.9k), which stays one leaf -> open pool precision 0.69, recall 0.857, binary accuracy 0.90, NMI 0.43-0.44
  (m8c39 0.73-0.75 / 0.71-0.75 / 0.90 / 0.39; ceiling 0.77 / 0.88 / 0.93 / 0.53). The r2 routing loss is gone; precision is limited
  by the root state's ~2.1k faint closed contaminants. M = MF on the open pool (the floor only merges closed leaves). Replicate
  Jaccard 0.992. arm_M.txt / arm_MF.txt.
- Arm M readouts (2026-10-24): merge-up units = open (the open leaf) + closed units, scoring equal to the leaves. About 1,700
  particles per seed are ceiling-reducible (1,296 of 2,131 contaminants, 409 of 792 lost open); the contaminants are faint (SNR
  0.0053 vs 0.0107). FINAL E-STEP (estep.sbatch: root stage-B project + leaf states -> one K-state reconstruct3D -> one refine3D pass
  with stage B options, m8c52t, 4.3 core-h per seed): precision 0.69 -> 0.74-0.75, recall 0.85-0.87, accuracy 0.92, NMI 0.49
  (consensus-pose ceiling 0.77 / 0.88 / 0.93 / 0.53). It fixes about half of the reducible contaminants.
  BUT the seeds agree less after it (open-pool Jaccard 0.873 vs 0.992; state NMI 0.30 vs 0.42). By direction: out-to-closed moves are
  right 85 % when both seeds agree (73 % one seed); into-open moves are coin flips (47 % / 24 %). The INTERSECTION of the two seeds'
  E-step open pools is best and reproducible: precision 0.780, recall 0.839, accuracy 0.926, NMI 0.504 (estep_agree.py).
  GENERALISES to the junk toys (estep_toy.sbatch on fulltoy_m52M_cd1 / m52MF_cd1, ~10 core-h each; scored by toy_junk.py +
  factorial's validate_tree3): ~1.4k truth junk move into the junk leaf, 245 / 241 out; net junk in genuine states -1.2k; junk
  precision 0.998, recall 0.44 -> 0.47. validate_tree3: MF NMI 0.284 -> 0.351, j2 0.21 -> 0.42, closed leaf 0.779 / 0.358 ->
  0.821 / 0.742, A-pure share 1.000 -> 0.923; M NMI 0.275 -> 0.339, j2 0.21 -> 0.45, closed leaf 0.776 / 0.358 -> 0.828 / 0.736
  (A mixed in at 0.77); the null-type state pairs mostly stay.
  It does NOT transfer to OLE real data (free K-state E-step erodes the reproducible root split): see [[final-estep-ole-real-data-2026-10-24]].
- STANDING (user, relayed 2026-10-23): "not sure why we would run cahrah with c1 and not c11 in the first place". CAHRA runs are
  C11 ONLY from now on; the C1 results stay as the record of why symmetry matters. The m8c26 C1 run (54132284) was cancelled; the
  m8c26 C11 run (54132283) is the current full-method test.
- Lesson: low-SNR symmetric heterogeneity needs the tree run with the point group. Run symmetric datasets in their point group
  from the start, not C1 first. C11 tree rerun: job 54119722 (tree_c11_s1/s2).
  Related: [[flexible-target-shortlist-empiar-2026-10-22]].
