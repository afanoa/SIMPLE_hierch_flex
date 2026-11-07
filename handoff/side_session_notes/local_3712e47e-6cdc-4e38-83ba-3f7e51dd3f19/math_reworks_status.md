# hetero_tree mathematical reworks: status (session tree-math-reworks, worktree tree-math-reworks-30ff30, branch claude/tree-math-reworks-30ff30 off claude/tree-dist 9005a8c77)

Work dir on batch2: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks (builds, gates/, dump_sites/, logs/).
Note on the branch: doc/implementation_notes/hetero_tree_reworks_2026_11.md.

## 2026-11-01 evening
- Read the audit (incl. the 6.7 amendment: junk screen stays, soft uses) and the site_rank / state_mergeup / controller code.
- Rework 1 (6.3) step 1: site_rank now dumps the hold-out-band Gram matrices (even.even, odd.odd, even.odd) of the endpoint
  deviation maps and of every final state's seed set (scored by posterior, no EM) to site_gram.txt; dev_gram added to
  state_mergeup. Decisions unchanged. Build build_dump (fast gate 15/15).
- Prototype: scratchpad subspace_proto.py reads site_gram.txt + the halves verdicts; rules: Fisher-z self reproducibility of the
  cross-half deviation, residual (principal-angle) test for collinearity, greedy complete-linkage grouping, seed-set tests.
- Running: 32 saved-site reruns (qtoy2_keep_s1, qtoynull_s2, olenull_s2 roots; the other trees had their masks pruned) and the
  five baseline gate trees (base_qtoy2 md1, base_r131 md1 seed 2, base_qtoynull/olenull/micnull md0) on build_dump with
  SIMPLE_TREE_KEEPMAPS=yes, keys = the recipe18d gate runs'.
- Sibling session node-mask-and-site-disco: no overlap (it swaps the INPUT half maps for 2D sub-classes); told me the r131
  fixture's cheap-children mask covers 93 % of the molecule support = the live failure case for 6.1.

## 2026-11-01 late evening
- Rework 1 (subspace site test) is built (build_sub, gate 15/15) and committed provisionally (aa42f6ccd, to be amended if the replay
  asks). Same-input rerun of 32 saved sites: nulls identical (0 poles both ways), qtoy2 sites 1-2 identical, sites 3-6 differ (more
  poles kept; tiny rest states of 30-88 seeds appear), 2.2-3.3x less wall, 4 EM fits instead of 28 per site. Gate trees: base arm
  (build_dump) and sub arm (build_sub) running; base nulls micnull / qtoynull pass (one leaf).
- Rework 2 (single-pass split test + numeric split_stats.txt + quality flag) is written on top (uncommitted), built (build_r2, gate
  15/15); its null trees are running (r2_*), the rest after the base arm frees cores.
- Prototypes 6.1 / 6.2 done (mask_disc_proto.py): 6.1 NO-GO (between-state variance is above the noise ceiling on 75-100 % of the
  sphere, null included), 6.2 partial (SVD yes, one share threshold stays). Details in the branch note.

## 2026-11-01 ~21:30 (batch2 time)
- Baseline arm done (nulls one leaf each; qtoy2 7 leaves NMI .712; r131 closed corner at r1_s2 .99/.51, as recipe18d OFF). Sub
  (rework 1) arm: nulls one leaf each, qtoy2 / r131 / ole finishing. Rework 2 nulls pass (one leaf) but that build skipped the
  production choose routing (SIMPLE_SPLIT_ROUTE default is forced to 'choose' in new()); routing ported into test_once, build_r2c
  compiling, the r2 arm is relaunched on it when ready. A scancel of the routing-less r2 qtoy2/r131/ole trees was denied by the
  permission classifier, so they run to completion (~75 core-h) and are reported as a no-routing variant.
- Bug fixed in test_once (weakest kept pair statistic overwritten by the threshold; merge-up margin only).

## 2026-11-01 ~22:15
- Rework 1 variant A: qtoy2 NMI .764 (base .712), nulls pass, r131 CORNER LOST at md1 - same two root poles as the baseline, but the
  seed-set z ranked site 3 over site 1 (closed population) -> different root seeding -> r1 re-test z 2.45 < 3.71. Variant B ranks by
  the weakest endpoint hold-out z (restores site 1); builds sub2 (rework 1B) and r2d (1B + rework 2 with choose routing) compiling,
  then sub2 r131 s2/s1 + qtoy2 and r2d qtoy2/r131/nulls. Seed-1 base/sub r131 replays also running.

## 2026-11-01 ~23:30 (batch2 time)
- Rework 1: variant A (seed-set z ranking) loses the r131 corner at md1 through the strongest_pole site order; variant B ranks by the
  weakest endpoint's hold-out z (sub2 arm running); variant C adds a 500-seed floor + fold of smaller states after a stage-A crash
  on a 4-particle rest state in the qtoy2 r1 node (sub3 arm running). qtoy2 variant A: NMI .764 vs base .712.
- Rework 2: with choose routing (r2c, on variant A sites) nulls pass, r131 r1 passes (corner delivered) and qtoy2 over-splits two
  combos by quality; cause = unequal N after the EM (old test: pairs at the small EM group, null at the stage-B N -> refuses; first
  single pass: big group at full N -> passes). Final form draws pairs AND null at the smallest EM group's N (build_r2f, replaying
  qtoy2 / r131 / nulls now). Seed-1 r131 replays of base and variant A also running; OLE scoring pending.

## 2026-11-02 ~00:30 (batch2 time)
- Rework 1 final form (variant C: subspace test + endpoint hold-out ranking + 500-seed floor) committed as 53593cf3f (amended). Variant B
  delivered the r131 corner at both seeds (base-identical tree at seed 2); variant C's own replay (sub3) and the one-N single pass (r2f,
  on top of C) are running; OLE arms for both are running too. OLE base vs variant A at md1: landscapes as close as the null.
- 6.4 prototype found one fold case (OLE r_s3/r_s5: Spearman -0.60, map angle 144 deg, 34 % reproducible residual): both tests agree,
  the map test adds information but needs an angle threshold. Rework 2 changes sit uncommitted on top of the rework 1 commit.

## 2026-11-02 ~01:30 (batch2 time)
- Rework 2 verdict: NO-GO as a default. The one-N single pass (r2f) passes the nulls but refuses r131's r1 at seed 2 (T 0.109 vs
  0.137) where the old test passed (0.154 vs 0.136) with the SAME EM and the SAME equal N; the null threshold does not depend on N;
  the verdict is set by the single random equal-N subsample draw. Remedy (not done): average T and the null over several draws.
  qtoy2: 9 leaves NMI .728 vs .707 (old test, same upstream). Seed-1 r131 of the one-N form running for a second data point.
- Rework 1: variant C (500-seed floor) passes r131 2/2 but changes small-node seeding on qtoy2 (NMI .707); variant D (B + fold below
  refine3D's 6 only) building/replaying as the intended final form (B's qtoy2 .764, r131 2/2, crash fixed).

## 2026-11-02 ~00:20 (batch2 time)
- Rework 1 = variant D committed (amended): nulls pass, qtoy2 NMI .764 (base .712), r131 corner 2/2 seeds (base-identical at seed 2),
  site_rank 2.2-3.3x faster. GO, pending the user's view on the ranking statistic.
- Rework 2 committed as implemented but NOT the default (callers on the production test): the one-N single pass loses r131 r1 at seed 2
  on both C and D upstreams by the single subsample draw. quality flag + split_stats + ptcl_weights written by state_split.
- Item 13: p_empty = 0 everywhere; the qtoynull root split is NOT a quality split by the truth (snrq d .03), so the flag correctly
  stays off; stage-B corr tracks snrq (r .7). Memory corrected.

## 2026-11-02 ~01:00 (batch2 time) - final state
- Branch claude/tree-math-reworks-30ff30: 6b4697c38 (rework 1, variant D, GO), ee764dda6 (rework 2, implemented, not default, NO-GO
  by the single-draw variance), dda3a18e5 + 3cbccf50c (note). Working tree clean. Nothing pushed.
- OLE md1: variants A, C and C+one-N give identical trees; all within the null of the baseline. Variant D's OLE tree pending (last
  data point); r2f seed 1 also lost the r131 corner (0/2 for the one-N pass vs 2/2 old test).
- Items 7 / 13 reported in the note sections 6.1-6.3; the quality flag stays off on qtoynull because that split is not a quality
  split by the truth (memory corrected). Campaign cost ~1350 core-h.

## 2026-11-02 ~03:00 (batch2 time) - draw-free statistic (overnight item)
- Branch rebased onto claude/tree-dist (cleanup e12be4e8d + variant D 26c350abb); my commits on top: 9c88430cd (rework 2, not default),
  two note commits, 9ee5d5c72 (draw-free pair statistic in the single pass; state_split runs it for the replay), note section 8.
- Draw-free T = per-shell cross-half product of full-N difference maps over (that + the groups' half-difference noise scaled to the
  common n), 3-shell smoothed, peak in the hold-out band; null = same at full N from within-strata shuffles. No subsample anywhere;
  merged groups = weighted means of their maps.
- Python check on saved stage-B maps: deterministic, monotone in n, equals the actual full-N cross-half FSC at full N.
- Running on build_df (instrumented: production test logged, draw-free single pass decides): 24 standalone split replays at 4
  subsample seeds x 6 nodes (old-T spread vs draw-free T) and the gate trees df_{qtoy2, r131 s2, r131 s1, qtoynull, olenull,
  micnull, ole}. Variant D's OLE tree (sub4_ole) still running.

## 2026-11-02 ~03:45 (batch2 time) - draw-variance table in (note 8.3, commit 339b681b4)
- Over 4 subsample seeds the production T spreads sd .02-.04 at the marginal nodes (r131 r1 corner node: .074 .137 .166 .097, 1 PASS
  in 4; variant-D r131 r1 seed 1: 0/4; qtoy2 r3 2/4, r5 1/4). The draw-free T is identical at every seed and sits at the old T's mean
  (.117 vs .119 at the corner node); its null is ~3x tighter (sd .01), so every tested node passes it at z 4.6-10.9, including the
  qtoynull root (which the locality guard then merges, as before).
- Gate trees on the draw-free single pass running (df arm, 7 trees, launched ~02:55); micnull done: root PASS z 82, guard merged, 1 leaf.
  Readouts (leaves, qtoy2 NMI, r131 corner both seeds, nulls, OLE landscape) and the go/no-go follow when they finish (~04:30-05:00).

## 2026-11-02 ~05:00 (batch2 time) - draw-free gate trees in (note 8.4; OLE pending)
- Nulls 1 leaf each. r131 closed corner retained at BOTH seeds (z 10.9 / 7.1) - the production test's seed-2 retention is 1 draw in 4.
- qtoy2: same 5-way root and r1/r2/r4 splits as variant D plus one split inside combo 28 at r3 (old T .165 z 3.25 < 3.38; draw-free
  .182 z 5.6): 9 leaves, NMI .742 (D .764, base .712). r131: one extra split per seed, a 0.93-closed-arm state beside the corner.
- Verdict: the statistic is right (no draw in the verdict, hard gates pass with margin) but the NMI dip means the gates do not clearly
  support a default change -> default stays on the production test (680bda2d5 reverts the caller; test_once stays one call away).
  The switch is the user's decision. Final build + fast gate of HEAD running (build_final, job 54836456).
- OLE md1 on the draw-free pass = variant D's tree to the particle (8 leaves, W1 0.000; every split z 8.5-28); within the baseline's
  null. Final state: branch HEAD 680bda2d5 (+ note commits), built in math_reworks/build_final, fast gate 15/15, nothing pushed.
