---
name: recipe12-lever-classes-2026-10-28
description: "recipe11 lever equivalence: disk cleanup, size-scaled threads, root site batch exact; shared-memory small nodes and 16k rows/part CHANGE results (site discovery 7 vs 4 sites, ARI 0.01-0.64) -> dropped; recipe12 definition"
metadata:
  type: project
---

2026-10-28, optimization session local_a012d0f3 (builds /scratch/cluster_scratch/atifao/recipe11, recipe12; plain-lever build renamed recipe11p).
- Exact (kept): disk cleanup incl. prune_partials path fix (qtoy2 byte-identical labels, 1.7 vs 2.9 GB); site/split/merge-up
  threads by node size (site 32 -> 10-15); rank all sites at once on nodes > 200k (same jobs, first splitting site in order);
  plain maps as nodes finish (scheduling only). Distributed site ranks are bitwise reproducible (r3222 twice, ARI 1.000).
- NOT equivalent (dropped, default off): shared memory on small nodes - same K* but labels ARI 0.64/0.13/0.01 at 7k/13k/21k, enters
  at site discovery (r3121: 7 vs 4 sites); gain small (-26 % core-h at 7k, no wall gain). 16k rows/part: changes part decomposition
  > 64k, untested -> back to 20k.
- recipe12 = recipe9 (has NFS done-file fix) + disk + site batch + threads + plain-as-finish + equivalence merge-up ON +
  strongest_pole default yes; tested on the subtree set only ([[test-on-subtrees-not-full-692k-2026-10-28]]).
**Why:** site discovery depends on part decomposition - a reproducibility sensitivity, not just a speed result.
**How to apply:** never re-propose shared-memory nodes or a different rows-per-part without an equivalence check at site discovery.
- 2026-10-28 ~18:30 recipe12b (recipe12 with recipe9's split/merge-up thread rules) == recipe9+strongest BITWISE on the r221
  subtree (7 nodes: same children, seeds, K*, sites; 18 stage-B maps byte-identical; r2211 pole margin 6.403 at 8 and 32 site
  threads). recipe12's divergence came only from split threads 8 -> 11 (state_split not thread-invariant). Savings: core-h -19 %
  (site ranking -65 %), disk final 13.1 -> 2.0 GB, peak -82 %; wall unchanged (queue noise). recipe13 = recipe12b + route A,
  recipe14 = recipe13 + choose default on.
