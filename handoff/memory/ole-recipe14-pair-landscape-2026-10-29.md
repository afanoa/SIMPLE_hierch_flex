---
name: ole-recipe14-pair-landscape-2026-10-29
description: recipe14 OLE seed pair reproduces its landscape (W1 .044, as strong .041, better than ref6~ref7 .055); 5-way root carries it, deep splits stochastic
metadata:
  type: project
---
recipe14 OLE trees, seeds 1 and 2: production keys, maxdepth 4, nodepar 4. 4.8 / 5.2 h wall, 218 / 223 core-h.
Runs are in beegfs hetero_tree_20261028/ole_r14/r14_s{1,2}; scored with scripts/hetero_tree_ole_repro/landscape_score.py.
- Landscape at the leaves, s1~s2: W1 .044, extrema .058/.021. strong s1~s2 is .041 and ref6~ref7 is .055.
- r14 vs strong is .058; r14 vs the old refs is .063-.075.
- r14 and strong reach the a3 high end at .58-.61; the old trees reach .53-.55.
- The 5-way root (route A pole seeds) carries most of the landscape, so the null shuffled within depth-1 groups is weak
  (.050-.053). The depth 2-4 splits do not reproduce by particle (ARI <= .33).
- Leaf maps s1~s2: coverage .80, share .88. This is the best pair so far; against ref6 the share is only .33-.35.

**Why:** this is the user's yardstick: whole-landscape and extrema reproducibility over two runs per build
([[landscape-must-be-reproducible-and-splits-data-derived]]).

**How to apply:** use r14 s1~s2 as the recipe14 floor. A depth-1-group null under-scales any build whose root splits more
than 2 ways. Doc: hetero_tree_ole_reproducibility_2026_10_28.md section 11.3.
