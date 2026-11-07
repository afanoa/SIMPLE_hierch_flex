---
name: flex-floor-random-fill-and-best-round-patch
description: refine3D_states flex floor delivers the LAST retry round and random-fills the rest (not 5%: 80% random); ~/SIMPLE_floorbest patch delivers the best round; not comparison runs 2026-10-01
metadata:
  type: project
---

`place_states_with_population_floor` (simple_flex_pca_model.f90, hael f1c3eb330 2026-09-12): rounds peel
under-floor clusters for good and re-place MORE clusters (K += deficit) on LESS mass, so qualifying counts fall;
on failure it keeps the 6 largest clusters of the LAST round and gives every other particle a uniformly random
label (no rationale in commit or doc/algorithms/heterogeneity_analysis/flex_pca.md). Surplus qualifying clusters
already get nearest-target assignment, so random is not forced. refine3D_states re-scores every sampled particle
against all states (prob_neigh_mode=state, ref index (state-1)*nspace+proj), so labels are only a start.

`not` (452,386 ptcls, nstates=6): floor 10% -> 112,691 random (25%); floor 5% -> 5 rounds, qualifying 4,4,3,3,1,
360,099 random (80%) though round 1 had 4 qualifying and ~13k in the rest. User: keep defaults, random is OK for now.

Patch: ~/SIMPLE_floorbest (batch2 worktree, detached f9bbf72ce + edit, own build, NOT committed) records the round
with most qualifying clusters (earliest on ties) and delivers it on failure; logs `delivering round=N`.
Comparison (both nparts=10): default `not/REFINE3D_STATES_OUTPUT_default_n10t8` job 53536843 (8_refine3D_states)
vs patched `REFINE3D_STATES_OUTPUT_floorbest_n10t6` job 53537152 (nthr=6 to fit the 70 free cores).
Earlier attempts 53531914 (beegfs outage), 53533881 (5%), 53536457 (176 cores) were cancelled.

**Comparison result (2026-10-01 15:42), same command, floor 10%:** default rerun: round 1 qualifying=2, 3 rounds, randomized
258,614 (57%; the first 10% run gave 112,691 -> the floor outcome is NOT reproducible run to run), states 100k/85k/74k/71k/61k/61k.
Patched: delivered round 1, randomized 0, states 265,116/129,559/31,859/11,379/11,052/3,421 (2 big clusters + 4 below floor).
