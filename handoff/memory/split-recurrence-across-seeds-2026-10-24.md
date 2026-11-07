---
name: split-recurrence-across-seeds-2026-10-24
description: "Seed recurrence of hetero_tree splits (bootstrap-style support) does NOT vet deep splits: AUC 0.975 on CAHRA top splits but 0.667 cdisc / 0.604 below the root; null floor 0.27-0.52; real deep splits land in different branches per seed"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-24T16:03:53.491Z
---

2026-10-24, for the coordinating session (user-approved analysis); worktree laughing-carson-4e8ca1, routes/toy55 `split_support.py`.
Output: rnasep_klem/factorial/validation/split_support.txt; write-up in factorial_site_layer_results_2026_10_21.md.
- Support of split X = best over the other seed's splits Y of the child-pairing mean Jaccard on the shared particles, times
  |S1 & S2| / |S1|.
- Pairs: p27cd1 vs cd2 (cdisc), and CAHRA M, MF and m8c39 s1 vs s2.
- Results: 19 truth-supported splits (median 0.932, min 0.164) vs 17 null (median 0.366, max 0.519). AUC 0.771 overall, 0.975 CAHRA,
  0.667 cdisc, 0.604 below the root.
- No cut works: 0.6 drops every null but also all seven deep cdisc real splits (closed / junk, 0.16-0.47).
- Causes: a null halving still matches the other seed's coarser split at 0.27-0.52. And seeds put real factor splits in different
  branches: p27cd1 finds closed under A-absent, cd2 under A-present.
- Merged units, by contrast, DO absorb null splits: the merge-up collapsed every genuine null pair in the m52 M / MF toys.
**Why:** a "keep a split only if it recurs across seeds" rule would delete real deep structure in factorial data.
**How to apply:** vet null-type splits by merge-up or composition, not by seed recurrence. Recurrence only confirms strong splits.
Related: [[tree-leaves-not-reproducible-across-runs-2026-10-21]], [[landscape-must-be-reproducible-and-splits-data-derived]],
[[factorial-site-layer-port-2026-10-23]].
