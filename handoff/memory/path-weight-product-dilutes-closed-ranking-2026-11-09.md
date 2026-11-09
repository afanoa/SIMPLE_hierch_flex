---
name: path-weight-product-dilutes-closed-ranking-2026-11-09
description: Bishop-Tipping path product of responsibilities (root -> node) on production's evidence adds closed purity only below r131 (r1312 +18-21 %, subtree r131 +18-43 %, Jaccard .13 -> .18) and removes it at r13 / r131 (-8 to -28 %): the closed signal sits at the r13 split
metadata:
  type: project
---

node_levers/nl_path.py, nl_path_seeds.py; handoff/node-levers/PATH_WEIGHTS.md (2026-11-09, zero compute). P_X = product of the
stage-B responsibilities along the lineage; at matched core sizes vs the one-level R_{X|parent}: r13 -181 / -295 / -235 closed
at the .9 / .95 / .99 sizes (the root factor excludes closed-rich particles of median R_{r1|r} .61); r131 -201 / -271 / -355 (r1
and r13 factors); r1312 +180 / +119 / +70 (x1.51 -> x1.78 at .9); r1311 and r132 "gain" only through leakage of the upstream
closed factors. Subtree seeds at r131: path cores purer on every run and cross-run Jaccard .132 -> .184 (shared 158 -> 352).
One-level closed AUROC by level: root .63, r1 .75, r13 .80, r131 .58.

**Why:** the coordinator asked whether the hierarchy's path weight beats the node's own responsibility for confident cores.

**How to apply:** do not rank cores by the uniform product; rank by the most closed-informative level (r13) or weight factors by
their information. Related: [[responsibility-selection-beats-hard-child-2026-11-09]], [[confident-cores-replicate-with-the-split-2026-11-09]].
