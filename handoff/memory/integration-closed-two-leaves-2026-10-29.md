---
name: integration-closed-two-leaves-2026-10-29
description: USER VERDICT 2026-10-29: the 692k r17d closed-lineage leaves r1311_s2, r1312, r1311_s1, r2322 are 4 DISTINCT closed conformers; keeping them apart is correct
metadata:
  type: project
---
USER VERDICT (2026-10-29, judging the four 10 Å maps by eye): "those 4 maps are distinct closed conformers".

692k integration tree (recipe17d, integration/f692_r17d_s1). The closed lineage ends in four leaves:
- r1311_s2: the corner, closed 1.00 / open 0.58.
- r1312_g1_2: 0.84 / 0.56. Its refused stage-B states were 0.72 and 0.94.
- r1311_s1: 0.77 / 0.57.
- r2322_g1_2: 0.78 / 0.69.
The merge-up (equivalence rule, mdelta 0.000075) never evaluated these pairs: their re-test stats exceeded the shuffle
threshold at every level, so they count as shown distinct.

From the halves and the 10 Å views (route9/integ_closed; pair_maps.py, view_cmp.py in split_merge_replay):
- Within a group (corner vs sibling; r1312 s02 vs s01) the difference is closed-arm occupancy only.
- corner vs r1312 carries a second difference a median 22 Å from the closed arm. It is present even in r1312's 0.94-closed
  state.
- r2322 differs about 30 Å from the arm.
- The read: r1312 is a distinct closed sub-state, not the corner diluted by open particles (its open-arm occupancy is
  unchanged).

**Why:** the user's eye confirms what the halves showed, so the tree's separation is right. Contrast the earlier r111-r1 pair,
which the user called close and "not worth a state" ([[close-closed-substates-not-worth-a-split-2026-10-28]]).
**How to apply:**
- Never fuse these four leaves; a merge rule that would fuse them is wrong.
- Count all four as retained closed conformers when scoring closed-state retention. The scorecard's "closed corner"
  (>= 0.95 / <= 0.6) admits only r1311_s2, so it undercounts.
