---
name: hetero-tree-flex-singleton-state-aborts-sta0
description: "TRAP: flex_pca's occupancy floor keeps a 1-particle runner-up state when fewer than 2 states clear it; the controller hands it to stA0, refine3D prob startup aborts (needs > PROB_MULTISTATE_MIN_POP = 5) and the node's particles end in no leaf"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-24T20:44:45.578Z
---

2026-10-24, RNase P junk-off tree full692_m33_s1 (job 54181559, m8c36). Both children of r22111, r221111 (4,801) and
r221112 (3,262), failed at stA0: "multi-state startup has an insufficient state population; state 2 has population 1".

- The chain:
  - flex finds one macro-cluster, then cuts 3 continuum k-means regions; two of them are single outliers (4,799 / 1 / 1).
  - The occupancy floor (simple_flex_pca_state_service.f90): when fewer than 2 states clear min_neff, it keeps the best
    state plus the runner-up "when it holds mass", so a neff-1 singleton survives as state 2.
  - case(S_FFLEX) passes every present label to stA0, and refine3D aborts.
- A failed node is not a leaf. finish() labels its rows 0 ("other" in hetero_tree_unplaced.txt), and they are in no merged unit.
- site_rank already drops pole seeds under PROB_MULTISTATE_MIN_POP + 1. Unguarded: fflex -> stA0, and the two-group fflex
  fallback that reuses the same labels.
- A third node, r112111 (5,203), failed the same way, for 13,266 rows so far.
  - 7 of the tree's 74 continuum-branch nodes (all 3.3-5.6k, depth 5-6) got a state of 5 or fewer particles.
  - Only those with TWO singletons aborted, because the floor kept a singleton as the runner-up.
- COMMITTED (user "fix the bug", relayed) as d3513ff70 on 16e52fbdd. Guard only; snapshot m8c61t = HEAD + guard,
  sources md5-identical, gated; r221111 replay places every row. A 4th m33 node (r1221121) had aborted: 19,322 unplaced.
- The guard was first built in snapshot m8c59t (= m8c58t + guard). Re-entry ([[tree-junk-pole-reentry-2026-10-24]]) stays uncommitted:
  - S_FFLEX unassigns flex states with <= PROB_MULTISTATE_MIN_POP members for stA0, and logs it;
  - the two-group fallback needs both groups above that minimum;
  - S_STAGEA and S_RECB fold small stage-A states into the largest (fold_small).
- Validated (3.0 core-h, stree_port/stree_reent_rnp2.sbatch): all three nodes run past stA0 and become one leaf each with every row placed.
  - The stage-A fold is untested on real data.

**Why:** a node lost this way silently drops its whole subtree's particles from the delivered labels.
**How to apply:** after any tree, grep tree.log for "FAILED (exit code" and the unplaced count before scoring it.
