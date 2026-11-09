---
name: textbook-soft-refit-r13-r1-no-clear-win-2026-11-09
description: textbook (M-step-only) soft refit at r1->r13 (5 seeds) and root->r1 (3 seeds): r13 same split as hard (+.03 concentration, latent .744 vs .733); r1 a different, more reproducible split (x1.92 vs 1.79-1.90) but lower recall and a worse latent (.65 vs .72); no clear win, 300 core-h
metadata:
  type: project
---

node_levers SOFT_LOOP.md 2026-11-09 (replays r13_b_i5*, r1_b_i5*; install5). r13: hard x1.70-1.71 / .354, soft x1.73-1.74 / .355,
Jaccard .92-.93 between arms, .90 across seeds. r1: hard x1.79-1.90 / .355-.389 (Jaccard .70-.84), soft x1.92-1.93 / .354-.358
(Jaccard .83-.84), between arms .17-.22; production r1 x2.12 / .380 remains best. Weighted fits show no shrinkage (corr .04-.07).

**Why:** the user wanted the soft refit tested where the parent cut loses closed rows (.68 kept at r1 -> r13).

**How to apply:** the soft refit is not a lever at r13 or r1 either; stop the thread unless a new mechanism is proposed.
Related: [[textbook-flex-weighting-production-level-2026-11-09]].
