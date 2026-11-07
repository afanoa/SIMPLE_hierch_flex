---
name: tree-pose-conditioned-retest-no-go-2026-11-02
description: "NO-GO: re-testing a hetero_tree split on its pose-stable particles (stage-B pose within 30 deg of the input pose) separates nothing: the draw-free T and the guard's L/L_ref of the stable subsets equal the full children's (ratio 0.99-1.01) on 14 nodes, real states, pose partitions (micnull, real PfCRT), quality and the 692k root alike; the children of a 'pose partition' are populations the EM found, the turn is a correlate; r4's difference from its siblings survives among stable particles (T .96, L/Lref .98-1.00)"
metadata:
  type: project
---
2026-11-02, session pose-conditioned re-test (branch claude/amazing-tesla-c2b47f on top of claude/tree-dist a9b6507a7; note
doc/implementation_notes/tree_pose_conditioned_retest_2026_11.md; scripts scripts/hetero_tree_proto_posecond/ with the result
tables in results/; data batch2 hetero_tree_20261028/posecond/cases/<case>/; 26.5 core-h, 22 jobs of 16-32 cores).
- Method: per node the children's full / stable (angle <= 30) / turned subsets reconstructed at the stage-B poses (reconstruct3D
  nstates=K box_crop 88 ml_reg=no, the tree's recB recipe) from a copy of stageB/b.simple with the labels written into the state
  slot (posecond_setstates.py; prg=selection REFUSES node-sized projects: STK/PTCL consistency check); then pair_stat_df at the
  node's logged n and pair_locality (Python ports), node mask and sphere mask. Full condition reproduces the tree (T .966 vs
  logged .955 etc.).
- Result (14 nodes, 30 pairs): T_stable/T_full 0.99-1.01 everywhere: mic_het root (real+pose) .970/.966, micnull (pose null)
  .876/.880, real PfCRT roots .975/.976 and .987/.989, depth-1 real .910/.908, quality .167/.165, r131 .503/.503, 692k r4 pairs
  .955-.963 vs .957-.964. Turned subsets carry the difference too (.89-.99). Within-child stable-vs-turned T .40-.79 < between.
  Guard L/L_ref unchanged across conditions (micnull .24/.24/.28). Mask irrelevant (same picture with the sphere).
- Why: stage B assigns particles by fit to the state maps; any subset of a child reproduces its map; the turn is a correlate of
  the worse-fitting population, not what defines the children. The micelle 'pose partitions' are reproducible population
  differences (micelle bins / quality), diffuse (micnull) or made local by the user mask (PfCRT).
**How to apply:** do not build a pose-conditioned gate; keep the pose numbers report-only (dflip, node flip fraction) and keep
the locality guard + effect floor as the only split gates. For the 692k r4: its difference from r1/r2/r3 is a diffuse population
difference (quality blur), not an alignment artefact; the stage-B-corr quality flag is the lever, not poses. The
posecond_setstates.py trick (float32 slot 23 of segments 3/6) is the way to relabel states in node projects from Python.
Related: [[tree-pose-artefact-test-2026-11-02]], [[micelle-hetero-fixture-mic-het-2026-11-02]], [[pose-flip-stats-in-split-stats-2026-11-02]].
