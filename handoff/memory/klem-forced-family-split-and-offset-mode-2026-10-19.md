---
name: klem-forced-family-split-and-offset-mode-2026-10-19
description: "why the klem tree overcounts states: NFAM=2 is a forced median cut (always ~24/25 cells) and per-family refinement makes the two families reproducibly different; plus the board selects on an image background/offset mode (uniform solvent drop, >25 A) that reproduces across independent trees. Root split test + board hp cut under test"
metadata:
  type: project
---

Null pair (one conformer), 2026-10-19: every motion-region board returns ONE state = its whole family, so the 12 leaves
are 2 family maps repeated; no particle attribute predicts the family (all AUC 0.49-0.51, logistic CV = majority).
Root split test (split_test.sbatch: parent cells -> 6 candidate groups, each reconstructed from the COMMON consensus
poses, cross-tree merge with noise-null threshold): two-conformer pair -> groups follow the conformers (c18-rich 0.67-0.89,
c49 0.93-0.98) but each conformer keeps one extra group; null keeps a ~15 % group at difference FSC 0.90. Its difference
map (diff_diag.py) = uniform negative solvent + higher particle contrast, power at shells 1-10 (>23 A) 100-3700x the
consensus, identical in both trees; no attribute separates it (AUC 0.49-0.51): selection on a per-particle background/
offset noise mode, deterministic across datasets. Merge-side background terms (NUIS=2) leave 3 groups at 0.26 vs 0.24.
Board fix: `hp <A>` option (patch_hp.py, klem2 build 12:50): shells coarser than hp get infinite noise in the commander,
the distributed master and the worker. parent_hp.sbatch reruns parent boards with hp 60.
**Why:** a tree must not invent splits: the split count per node has to be data-derived (user: real landscapes are not
balanced trees). **How to apply:** per node propose <= 6 groups, test them on common poses across the two trees, keep
K* (1 = leaf); never refine children before the split is accepted. Related: [[klem-merge-threshold-needs-noise-null-2026-10-19]].
