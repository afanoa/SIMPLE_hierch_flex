---
name: ole-timing-triple-crops-verdict-2026-10-27
description: "OLE tree A/B (current / crops / crops+sitebatch+endphase+threads): crops 0.68x wall and core-h per node, root split kept; leaves differ only via r1's run-to-run split flip"
metadata:
  node_type: memory
  type: project
  originSessionId: 0f67cd42-16a3-4392-a859-a59fdac654fc
  modified: 2026-10-28T02:37:18.773Z
---

The OLE timing triple ran on 2026-10-27: trees ole_nosl2_{base,opt,full} under /mnt/beegfs/elmlund/afan/tree_port_optim_20261027/trees (seed 1, nodepar 8, site_layer=no, ITA = ITB = 1).

**Crops** (fflex delivered maps at BOXA, seeds/recB cropped, postprocess=no where only halves are read):
- Wall 557 -> 378 min.
- Core-h per node 13.1 -> 8.9 (0.68x).
- Root partition kept (ARI 0.97).

**The "full" arm** (crops plus sitebatch, endphase and threads) took 474 min.
- Its longest chain matched opt's through depth 5 (321 vs 321 min); its tree then split once more.
- The merge-up fell from 21 to 9 min.

**Leaves:**
- Leaf ARI between arms was only 0.25-0.27.
- r2's subtree reproduced (depth-2 ARI 0.91-0.97). r1's split flipped in every pair, even between opt and full, which had identical inputs.

**Why:** a tree A/B on OLE cannot use leaf identity as its gate; r1 is a run-to-run site-order pick (see [[ole-seed-split-choice-2026-10-25]]). Compare at the root and r2, or use a same-build repeat for the floor.

**How to apply:**
- Quote per-node core-h from sacct, not step times: queue waits (median 14-120 s per job) inflate step times, and they made "plain" look 4x slower.
- The end-phase tail is now plain (queue-bound) plus final_maps' 4 reconstructions run one after another (about 15 min).
- Related: [[m33-tree-wall-profile-2026-10-25]], [[tree-runtime-goal-5-7h-at-700k]].
