---
name: tree-latent-verdict-2026-11-07
description: "VERDICT: one shared sparse-PPCA site basis + a tree-structured XD mixture + sat5 fission splits beats the hetero tree on qtoyflex, ties it on qtoy2, gives 6x at IgM 10 % (.26 vs .042) but < lever 3 (.338); 3-state arms need K=3 proposals"
metadata:
  node_type: memory
  type: project
  originSessionId: 1218bb2b-a21a-4d25-9e48-c40ddd89ceff
  modified: 2026-11-08T02:02:25.638Z
---

The tree-latent session (2026-11-07; batch2 `handoff/tree-latent/VERDICT.md`, runs `math_reworks/tree_latent/`) built the cohesive
tree in latent space. It shares one site basis tree-wide, deconvolves per particle (posterior precision A), proposes splits by
K=2 XD, and accepts a split only on a sat5 DISCRETE call.

| Fixture | Latent tree | Comparators |
| --- | --- | --- |
| qtoyflex | 3 leaves (arm whole, holo, closed); closed F1 .558, holo .713 | hetero tree .39-.47 / .62 |
| qtoy2 | NMI .75-.76; closed and holo isolated | pp_on .811; the conformer 1 \| 2 pair merged |
| IgM 10 % fix | five-arm product: joint .260, forest .258 | hetero tree .042, lever 3 .338 |
| IgM 5 % fix | one leaf (z_mid < 3) | |
| Nulls | one leaf everywhere | |

**Why:** the user wants "a cohesive model over the entire tree". This prototype says the model works and that the K=2 split rule
limits it. A balanced K=2 split of a 3-state arm leaves a two-state side that sat5 reads as a continuum.

**How to apply:** the next step is K=3 proposals with a pairwise saturation test, in the forest form (one tree per site). Keep
memberships on y1 ([[fission-select-members-on-y1-labels-may-read-y]]). Related: [[data-fission-merge-test-verdict-2026-11-07]].
