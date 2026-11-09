---
name: closed-is-f1-factor-across-root-branches-2026-11-09
description: "VERDICT: RNase P closed is an F1 (lower region, (30,-43,-41) A) factor crossing the root split: lf and r17d roots cut along F2, the closed-concentrating depth 1-3 splits cut along F1 in r1, r2, r5 (r1/r2 same shift cos .88-.99, r5 another direction); site masks all pruned in both trees, leaf maps used. F1 x F2 two-way product: states factorial (free table z <= 2.1), ARI .28 vs .23, F1 closed side only 5/7 seeds: pre-set bar not met"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f195b76-225b-466c-a432-7fc8360148e8
  modified: 2026-11-09T19:26:36.619Z
---

Coordinator task (user: "Yes, run both"), 2026-11-09.
- Notes: `hetero_tree_20261028/handoff/ole-panel/CROSS_BRANCH_F1F2.md`.
- Scripts: `ole_panel/scripts/cb_check2.py`, `cb_leafmaps.py`.
- Outputs: `ole_panel/rootsp2/cb/`.
- Cost: 6.2 core-h of analysis, no SIMPLE runs.

**What survives.** Every site mask is pruned in lf and r17d: 463 + 434 entries, including `fsite_*/refine_mask.mrc`. Interior nodes have no maps. The leaf maps (with halves), each node's `rows.txt` and `sites.txt`, and the factorial pole maps survive.

**Check 1: closed crosses the root split.**
- The root splits on F2 (r1 and r2 at AUROC .68-.81).
- The closed-concentrating splits below run on F1, the lower region (x2.3-7.0, AUROC .64-.78): r1->r13, r21->r211, r211->r2112, r5->r51, r511->r5111; in r17d, r1->r13 and r13->r131.
- Closed leaves of r1 and r2 shift the same way on F1 (cos .88-.99). r5's shift another way (cos -.16 to .62).
- The leaf-map contrasts sit on F1 (cos .43-.64) and not on F2, F4 or the r131 block (<= .21). r1 vs r2 contrast cosine .85.

**Check 2: F1 x F2 two-way product.**
- The states are factorial: the free 2x2 table never beats the product by more than 2.1 SE.
- The joint 4-D mixture beats the product by z 6.6-9.9. Read this as possible shared noise: only marginal SDs exist, so cross-block noise is not deconvolved.
- Product ARI .28 against today's root split .23 (minimum .06 against .13).
- The F1 closed side carries closed in only 5/7 seeds (x1.23-1.40); in seeds 5 and 6 F1's bimodal axis is not closed.
- Best product cell x1.29-1.85.

**Why:** this tests the K-way-root-by-sites redesign.

**How to apply:**
- Promote F1 by where closed lives, but a plain root-level two-way F1 split is unreliable. Try a K-way or F2-conditioned F1 call.
- Masks for lineage checks must be kept at run time; production prunes them.
- Related: [[root-sparse-sites-seed-test-2026-11-09]], [[site-coupling-factorial-verdict-2026-11-07]].
