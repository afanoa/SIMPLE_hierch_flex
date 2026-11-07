---
name: hetero-tree-math-reworks-setup-2026-11-01
description: where the math-rework builds, gate replays and site reruns live on batch2, and the replay recipe (fixtures, keys, masks-pruned trap)
metadata:
  type: project
---
2026-11-01, session tree-math-reworks (worktree tree-math-reworks-30ff30, branch claude/tree-math-reworks-30ff30 off claude/tree-dist
9005a8c77). Everything on batch2 under /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/: builds build_dump (= 9005a8c77 +
Gram dump), build_sub (rework 1, subspace site test), build_r2 (+ rework 2 single-pass split test); gates/<arm>_<ds> replay trees
(arms base / sub / r2; ds qtoy2 md1, r131 md1 seed 2 from ancestral/fixtures/r131, qtoynull / olenull / micnull md0, ole md1),
launched by launch_gates.sh <install> <arm> [ds..] (keys = the recipe18d gate runs', SIMPLE_SPLIT_CONTRAST=maps, KEEPMAPS=yes);
readout.sbatch ARM=<arm> writes gates/<arm>_readout.txt (null leaf shares, qtoy_leaves.py, closed_nodes2.py). dump_sites/ and sub_sites/
= the same 32 saved sites re-ranked (site_rerun.sbatch over a list of "<src site dir> <out dir>").
TRAP: only trees run with SIMPLE_TREE_KEEPMAPS=yes keep sites/site_*.mrc; r11_new9, r131_OFF_s2, micnull_PG_s1 had theirs pruned, so
site_rank replays there need a new tree. qtoy2_keep_s1, route9/nulls/{qtoynull,olenull}_s2 keep them.
**How to apply:** reuse the arms and launcher for any further tree change; one site_rank at 60-140k rows takes 20-120 s on 16 threads.
Related: [[hetero-tree-algorithmic-audit-2026-11-01]], [[hetero-tree-prune-deletes-split-inputs]].
