---
name: soft-carry-no-gain-on-fixtures-2026-11-09
description: carrying refinement responsibilities down the hard tree (Bishop-Tipping path weights) gains nothing on truth or reproducibility on qtoy2 (ten trees) and qtoyflex (three seeds); lost states leave confidently; completing missing factors trades purity for recall at a loss
metadata:
  type: project
---

node_levers/nl_softtree.py, handoff/node-levers/SOFT_CARRY_FIXTURES.md (2026-11-09). Strict path arg-max = hard leaves (> 99 %
on qtoy2, exactly on qtoyflex's whole-node leaves); expected soft recall .3-1.2 % below hard on qtoy2; the completed variant
(missing factor = child's row share) moves 8-9 % of qtoyflex particles, combo 18 recall +.07-.09 but purity -.12 (F1 -.07 to
-.09), NMI -.01, cross-seed ARI -.02. Losses at cuts are confident (2 : 1, R p50 .02-.2) except qtoyflex r1 -> r12 for combo 18
(R p50 .76). Sharpness: qtoyflex deep nodes .81-.88 median with 21-30 % < .7 (softer than RNase P's .86-.92); nulls .81-.83.
IgM skipped (stock tree captures no arms at the root; 520-1,350 core-h per tree).

**Why:** the user asked whether soft carrying beats hard where states overlap; the hard tree holds no evidence off-path, so a
carry can only re-weight within a particle's own branch.

**How to apply:** do not pursue soft label carrying on hard trees; the responsibilities' value is ranking / selection
([[responsibility-selection-beats-hard-child-2026-11-09]]). Related: [[path-weight-product-dilutes-closed-ranking-2026-11-09]].
