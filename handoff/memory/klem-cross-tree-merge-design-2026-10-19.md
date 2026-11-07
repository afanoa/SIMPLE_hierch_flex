---
name: klem-cross-tree-merge-design-2026-10-19
description: within-tree half maps cannot judge state differences (selection bias is in both halves); the merge/support test must compare two INDEPENDENT trees (disjoint particle halves A/B): merge_cross_tree.py + cross_tree_test.sbatch; synthetic validation running
metadata:
  type: project
---

Finding (null 200k, 2026-10-19): the 15 same-conformer leaves survive view matching (14-15), per-shell spectral
matching (13) and rigid-body docking (15; every leaf docks with the same 0.9 deg / 1 px offset to the consensus).
Cause: a board selects a leaf's particles by their match to its own noisy cell map (cross-half E-step converges to a
fixed point where both halves share the selected noise), so the even AND odd maps of a leaf carry the same bias and
the difference between two leaves reproduces between halves without being conformational. No post-hoc fix on one
tree's half maps can detect this.
Design: gold-standard trees. Split the particles into disjoint halves A and B (split_halves.py: alternate rows ->
sim<MODE>_A/_B with marker files so run_toy55_scale skips simulation; poses refined independently), run the full chain
on each (synth_tree_prep_full MODE=<m>_A/_B), then merge_cross_tree.py RUN_A RUN_B cons OUT: Hungarian match of the
leaves by deviation correlation (15 A, union of region masks), unpaired leaves dropped, pairs p,q merged while
FSC(A_p - A_q, B_p - B_q) < 0.143 (best 3-shell band 40 A .. pooled cross-tree resolution; union mask or envelope);
writes OUT_A.npz / OUT_B.npz (labels; 0 = dropped) for consolidate_states.py on each half (lab 0 now skipped);
cross_tree_test.sbatch MODE= runs merge + consolidation + truth purity per half. Expected: null -> 1 state, two_50k ->
2, 55_200k -> fewer states than the 15/24/71 within-tree counts. Watcher chains null_200k, two_50k, 55_200k halves.
Cost: two trees of half size = about the same compute as one tree; real data (f692e) would need two 346k-particle trees.

**Why:** the user's question "are we only going to split as much as the data supports" has to be answered by a test
whose noise is independent between the two things compared; within-tree even/odd is not.
**How to apply:** every support/merge decision in the chain should become cross-tree; see
[[klem-null-leaves-sort-by-alignment-quality-2026-10-19]], [[klem-scaling-synthetic-55-2026-10-18]].
