---
name: mergeup-draw-free-port-2026-11-02
description: "the equivalence merge-up ported onto the draw-free single-pass test (c569be829, branch claude/optimistic-chandrasekhar-10506d = tree-lever-floor + port; build_mu on batch2): 692k r17d merge-up re-run 65 units vs saved 63 (the two borderline merges 43+47, 19+23 now distinct), 64.7 vs 73 min; six gates identical to bc; TRAP - trees older than c7ac3b0bd have no split_stats.txt, so the merge rule reads margin 0 there unless a shadow tree is given through SIMPLE_MERGE_TREEDIR (integration/shadow_split_stats.py)"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e9dba0b-f95e-43c4-bf4b-16ee45f9990d
  modified: 2026-11-02T20:35:00.000Z
---

2026-11-02, session tree-mergeup-port (worktree optimistic-chandrasekhar-10506d; the branch claude/tree-lever-floor-bccc29 lives in
its own worktree, so this worktree's branch was reset to its tip c5176bddc and the port committed there as c569be829; fast-forward
tree-lever-floor-bccc29 to it when accepted). Note: hetero_tree_reworks_2026_11.md section 11.
- The port: state_mergeup%run calls test_once with self%mtests (the tree's largest split test count -> z_crit); ugrp(unit) = the
  final index of the EM group the unit started (decide_once's gof through rank and the choose renumbering rkc), 0 when the EM
  drained it (passes alone); decide_once%shown_equal reads r_eq at the statistic's common n (= split_effect's n), not the pair's
  smallest member count. Deleted test/decide/eqn_stats/null_upto/pair_stat/blom_emax/rnd_subset + nmed/nmax/nnull, nul/jsc/l_df,
  MINSEL_JUNK/JUNK_CORR: 750 lines (5,325 -> 4,608). The commander is unchanged.
- VERDICT on the 692k (integration/mu_mergeup_r17d_s1, job 54866127, 48 cores): 64 min 44 s vs 73 min, 51.9 vs 58.5 core-h, 65
  final units vs 63. Both runs merge r2311_g1_2+r2312_g1_2 and r2122_s2+r2511_s1; the saved run's other two merges
  (r3112_s2+r3111_s1, r2211_g1_2+r2212_g1_2) passed the margin by 0.000006 / 0.000001 only at the pair's member count n under the
  old max-of-20-shuffles band (0.117); on the draw-free null (sd 0.008 vs 0.023, threshold 0.041) they are distinct at z >= 3.7.
  The four closed-lineage leaves are in no merge in either run. Gates: nulls 1 leaf, qtoy2 8 -> 8 NMI .812, r131 both seeds
  4 -> 4 with the closed corner, all identical to bc; OLE 8 -> 8 with 0 merges, landscape W1 0.000 vs bc, merge-up 17 min in both arms. Note commit on this branch follows c569be829.
- batch2: math_reworks/src_mu + build_mu (fast gate 15/15), gates/mu_* (bc = the reference arm: same tree code, production
  merge-up), integration/mu_mergeup.sbatch (TREE, BUILD, OUT, SHADOW=yes) = the merge-up alone on a finished tree in its own
  directory: labels back on a copy of mergeup/final.simple with `simple_exec prg=selection infile=hetero_tree_labels.txt
  oritype=ptcl3D` (7 s), then state_mergeup with the saved merge_tree.txt / leaf_masks.txt.
**TRAP:** since fe02cffb8 the merge rule reads only node_*/split/split_stats.txt; the Oct-29 r17d 692k tree has none -> margin 0,
no cannot-link pairs, nothing merges. shadow_split_stats.py rebuilds them from the production test's split logs with the deleted
node_effect's rules (59 nodes, margin 0.000075, 21 pairs = the 73-min run's inputs) into OUT/shadow_tree for SIMPLE_MERGE_TREEDIR.
Trees from build_lf on (f692_lf_s1) carry the numbers natively.
Related: [[tree-692k-wall-profile-2026-11-02]], [[equivalence-mergeup-verdict-2026-10-27]], [[joint-gate-c7ac3b0bd-2026-11-02]],
[[split-test-single-draw-variance-2026-11-02]], [[integration-closed-two-leaves-2026-10-29]].
