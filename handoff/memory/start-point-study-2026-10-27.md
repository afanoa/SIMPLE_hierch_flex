---
name: start-point-study-2026-10-27
description: "hetero_tree node start-point study (2026-10-27): state_anneal and majority-vote consensus FAIL; picking the most confident greedy restart (stage-B mean max state posterior) picks the truth-best basin; the qtoy-matched NULL splits by quality and state_split accepts it"
metadata:
  node_type: memory
  type: project
  originSessionId: 02205d1d-514c-444f-b666-b83e6423575a
  modified: 2026-10-28T01:39:26.089Z
---

2026-10-27, session start-point-sensitivity-cb0d4a, standalone refine3D fits of the node recipe (A3+B3), cluster folder
/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/start_point_v1. Write-up:
doc/implementation_notes/planned/hetero_tree_start_point_sensitivity_2026_10_27.md in that worktree.

- No existing knob softens the STATE under prob_neigh geom: prob_athres only shapes the projection draw; snhc is global.
- New opt-in key state_anneal (uncommitted): fails. r1211 latent start stays in its basin, site-pole start loses 7-18 pts.
- Not non-convergence: B9 runs stop converged after 4 its with unchanged partitions (distinct fixed points).
- Majority-vote consensus and cross-restart reproducibility follow the restart mix (r1211: wrong basin).
- Confidence selection (highest mean max posterior from prob_state_evidence.txt of the STAGE-B partition, greedy
  restarts) picks the truth-best basin: rank corr r1211 0.88, CAHRA 0.93, qtoy2 0.89. Fails at stage A on r1211.
- qtoy2 and CAHRA barely start-sensitive; RNase P r1211/r12111 and OLE r1111/r111 are.
- NULL (simqtoynull, one conformer): 7/8 restarts give the same snrq split (1.25 vs 0.82), confidence 0.86-0.90,
  halves gate 0.97; published state_split accepts it (pair stat 0.66 vs threshold 0.096, K*=2). Stratified null is
  view x defocus only. OLE-matched null not run (cancelled at wind-down). Study closed; 545 core-h.
- Iteration budget (same night, user approved ITA=ITB=3 for the integration build): 1+1 scatters the r1211 closed state
  over 3 leaves (best .691/.334); 3+3 keeps it in one leaf (.474/.914) but cannot isolate it without the closed-pole seed
  (reference tree: 3 site poles at r1211). qtoy2 3+3 depth-1 states = baseline spread (NMI conformer .508 vs .513).
- Below maxdepth, children get split-step regrouped groups (13-20 % moved, follow junk), not the node's stage-B states.

**Why:** the user wants reproducible landscapes; this shows what fixes the start and what does not.
**How to apply:** don't rerun annealing or consensus arms; restart selection only after the split test rejects the
quality null. Related: [[landscape-must-be-reproducible-and-splits-data-derived]], [[tree-design-feature-freeze-2026-10-25]],
[[gate-fakes-always-find-structure-2026-10-25]], [[hetero-tree-node-refinement-is-local-not-frozen]].
