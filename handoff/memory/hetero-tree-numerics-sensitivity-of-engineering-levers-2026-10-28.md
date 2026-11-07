---
name: hetero-tree-numerics-sensitivity-of-engineering-levers-2026-10-28
description: "Which hetero_tree speed levers keep RNase P results identical: distributed steps are bitwise reproducible; site_rank is thread-invariant; state_split thread count and shared-memory/part decomposition are NOT (1 row moved flips the next node's strongest pole)"
metadata:
  node_type: memory
  type: project
  originSessionId: c8c6f555-6cba-44a4-9ec4-48197b959f2f
  modified: 2026-10-29T01:21:35.093Z
---

Measured 2026-10-28 (tree-optim-budget session) on RNase P subtree r221 (recipe9 tree, 26.9k rows, 7 nodes), same inputs and seed:
- A distributed node replayed twice on the same build: identical splits (ARI 1.000). Run-to-run noise is zero at fixed settings.
- site_rank at 8 or 11 threads vs 32: identical seeds on every site, the same pole margin (6.403). Thread sizing of site ranks is safe;
  it cut site core-h 65 % (site ranking was 30 % of a 692k tree's core-h).
- state_split at 11 vs 8 threads moved 1 row of 26,900. The child node, starting from that one row, flipped its strongest-pole
  decision: margin 6.40 kept became 1.67 rejected, and the split diverged.
- Shared memory (nparts=1) vs 4 distributed parts changed site discovery (7 vs 4 sites); node labels at ARI 0.64 / 0.13 / 0.01.
- Equivalent levers (recipe12b, bitwise on all 7 nodes and on the r1211, qtoy null and OLE null gates):
  rolling disk cleanup; rank all sites at once; site threads by size; job CPUs = threads; plain maps as nodes finish.
- Row-local node projects under site_layer=yes (recipe15) were also equivalent on r111 route A at maxdepth 1: identical
  labels; flex at 1e-6 moved site masks by 1-6 voxels with no decision changed.
- qtoy2 did NOT show the split-thread sensitivity (labels identical). Equivalence checks must use RNase P nodes.

**Why:** the user's rule "we cannot sacrifice science for speed"; RNase P amplifies any float change at near-threshold decisions.
**How to apply:** keep state_split / merge-up / part-decomposition settings fixed unless gated. Check an engineering lever bitwise on
an RNase P subtree from copied parent files (recipe12/subtree.sbatch on beegfs). Related: [[numerically-close-not-bitwise]],
[[ole-repro-seeding-near-ties-2026-10-28]].
