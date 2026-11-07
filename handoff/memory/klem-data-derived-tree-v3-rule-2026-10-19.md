---
name: klem-data-derived-tree-v3-rule-2026-10-19
description: "data-derived klem tree that counts right on synthetic sets: per node NSEED board starts (keep highest mean max-responsibility), <=6 candidate groups, splitbench3.py (in-memory equal-N Fourier-insertion maps, coarse-band pairing, window stat, max-of-shuffles threshold, per-shell blend archetype pruning, soft-EM purification, re-test). null 1, near-null 1, two 2, two clades 3"
metadata:
  type: project
---

toy55 scripts (batch2 $W): dtree_v5.sbatch (driver; MODE, TAG, CONFS=a:b:.. root cut, NSEED, PARNODE, NMIN, MINGRP,
BENCHJOB=splitbench_node3.sbatch), dtree_board2.sbatch (SEED -> random start), split_candidates.py (MINSIZE), splitbench3.py
(cache benchcache_<MODE>_<H>_c48.npy), dtree_plot.py (tree figure), truth_dendrogram*.py (reference hierarchy of the 55).
Why each piece (2026-10-19): SIMPLE reconstructions gave shuffled-pair null up to 0.45 (bench: 99 % 0.07-0.09); one board
start failed in one half (clades not separated) -> keep the most confident of 3-4 starts (the board's own reproducibility
score picks the WRONG start); hard partitions imprint a filtered copy of the real difference on intermediate groups
(per-shell blend residual 0.45 -> 0.12, soft EM merges the rest); unpaired groups follow their nearest paired group.
Results (v3 at the root): null_50k 1, 52+53 1, single conformers 1, two_50k 2 (c49 0.94-0.95 / c18 0.82-0.83), orange+red
clades 3 (orange 0.90-0.95, red {5,6,21..} 0.98, red {47,43,7}+blend 0.67) = reference branches.
**How to apply:** use v3 for node decisions; check every new rule against null + shuffled controls first.
Related: [[klem-split-test-statistic-traps-2026-10-19]], [[klem-forced-family-split-and-offset-mode-2026-10-19]].

TRAPS fixed 15:00-15:40: (1) bench logs print "K* =" twice (board groups, soft EM): the driver must take the LAST;
(2) a trailing "# comment" on a bash line silently disabled the child-label copies that followed it; (3) child ids
were ordered by size per half, so node r.1 was green-rich in half A and purple in half B and every depth-1 test
compared unrelated particles (55-tree: 4 leaves, nothing split) -> children numbered by COMBINED size (splitbench3.py);
(4) equal N = smallest candidate group (582) kills power: MINGRP 1000-1500. Results after (1)-(2): two_50k tree = root ->
2 leaves (c49 0.94-0.95, c18 0.82-0.84); orange+red cut = 4 leaves (orange core, red {5,6,21} 0.98, red {47,43,7}
0.88-0.90, boundary {41,29,26}+{47,43}); null, 52+53, single conformers = 1.
