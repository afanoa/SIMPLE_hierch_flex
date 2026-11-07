---
name: lever3-latent-refs-seeding-verdict-2026-11-07
description: "VERDICT 2026-11-07 (IgM 10 %/5 % fixtures, true arm masks): latent-built references + arm-masked state_eval as a node seeding rule finds the arms far beyond the tree (5-arm product .256 vs tree .042 at 10 %), but is discrete-specific only for 2-state arms at round 1; 3-state splits and later rounds carve the null's continuous angles; Part 3 on sparse-PPCA supports is the best data-only chain (10 % product .338); Part 4: GAPw accepts the site splits (9/10 power, 0/10 false accepts) but fails on qtoyflex's curved continuum and closed state"
metadata:
  node_type: memory
  type: project
  originSessionId: 77c8a7fa-5cd0-4244-b81c-f308ec6cd639
  modified: 2026-11-07T22:02:50.799Z
---

Handoff: `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/igm-refs/VERDICT.md` (Part 1 true masks; Part 2 zh1 data-only masks). Runs: `math_reworks/latent/igm_refs/<tag>_a<k>`. Scripts: `igr_*`.
Rule: seed = T1's 3-D GMM on the arm-focused flex coordinate; references = the root's stage-B call with the seed labels forced; then
rounds of free state_eval (maxits 1) with the arm mask, each rebuilding its maps from its new labels.

- **10 %.** Product vs 72 combos: seed .222 -> R1 .235 -> R3 .256 (tree .042). Per arm: a1 .300 -> .351 (R1) -> .370; a3 .292 ->
  .346; a4 (failed seed .003) -> .139 by R3.
- **5 %.** Product .054 -> .098 (tree 1 leaf). a1 .000 -> .137 by R3.
- **Round 1, null.** 2-state arms gain <= .004 (pass). 3-state arms gain only on the null (+.015 to +.030), while the fixture's
  gain <= .02. Later rounds carve the null at any K (10 % null a1 .033 by R3). Null products: 10 % .059 -> .091, 5 % .021 -> .035.

- **Part 2 (data-only zh1 masks, done 14:09).**
  - 10 %: the product reaches .165 -> .225 (R1) -> .251 (R3), equal to the true masks (.256). a1 and a3 match the true mask;
    a4 beats it (R1 .085 vs .013, R3 .194 vs .139). The 3-state arms start weaker (seeds .11 vs .22) and recover to .17-.18.
  - 5 %: every zh1 seed is <= .001. The rounds recover a1 like the true mask (.131 vs .137 at R3), little else; product .044 vs .098.
  - The null behaves as with true masks: 10 % product .073 -> .093; 2-state null arms <= .007 at R1.

- **Part 3 (sparse-PPCA supports + their own block coordinates, fully data-only, done 15:45): the best chain.**
  - Product at R1: 10 % .338 (seed .333; true .235, zh1 .225, tree .042); 5 % .104 (seed .050; true .075, zh1 .025).
  - At 10 % the masked round nudges the 2-state arms toward its fixed point (a1 .394 -> .378) and lifts the 3-state ones (a5
    .295 -> .318). At 5 % it rescues arms: a3 .012 -> .139, a4 .011 -> .117.
  - Null: the sparse seeds already split the continuous arms (~.15 vs angle bins). Rounds add <= .007 on 2-state blocks at R1 and
    +.012-.034 on 3-state blocks.
  - Trick: resample the box-70 support centre-aligned to the 200^3 mask grid, core at .5, then T1's edge (dilate 3 A, Gaussian 6 A).

- **Part 4 (discreteness test, FINAL 17:20):** GAPw on the round-1 K=2 groups: 9/10 power, 0/10 false accepts (IgM arms 9/9,
  null blocks and the non-arm block rejected); fails on qtoyflex both ways. Details: [[gapw-site-acceptance-test-verdict-2026-11-07]].

**Why:** the coordinator's lever 3 asked whether reference-seeded, arm-masked assignment could seed tree nodes.
**How to apply:** if it is ported, use one round, 2-state splits, plus a map-level split test (GAPw / cross-half). Do not iterate
blindly: 3-state splits and extra rounds manufacture states from continua. Trap: the 3-D GMM seed is recipe-fragile; report the
seed NMI with every number. Related: [[igm-arm-flex-prior-verdict-2026-11-07]], [[stageb-oracle-references-limit-2026-11-07]],
[[igm-focused-flex-mask-is-bottleneck-2026-11-07]].
