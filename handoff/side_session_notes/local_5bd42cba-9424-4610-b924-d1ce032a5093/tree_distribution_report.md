# hetero_tree distribution: overnight report (session funny-pike, 2026-11-01)

Work dir on the cluster: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/dist_sdisc
Source: worktree branch claude/tree-dist at 21fac34db (fork hetero-tree) + uncommitted edits.

## Diagnosis (NOT tree 54765571, sacct + tree.log)
- 10.6 h of worker activity, 663 core-h, mean 62 cores.
- Every idle gap > 2 min (5 gaps, 2.25 h) is the master running SITE_DISCOVER in-process
  (r2 56 min, r21 45 min, r1 16 min, r23 10 min, r 8 min). The master is single-threaded and
  blocks every other node while it runs, so sibling nodes stall too.
- All other master work ("controller ... s") totals about 31 min over the run.
- nodepar was not the binding limit: the tree's own serial chain and the discovery stalls were.

## Changes (uncommitted, in the worktree)
1. Site discovery runs as a worker job (new step S_SDISC, `prg=site_discover`, 64 GB, nthr=nthr)
   in node/sites; the master only reads sites.txt when it finishes. Restart: S_SDISC re-runs from
   recB0. UI: site_discover accepts pgrp.
2. discover_sites, exact speedups: candidate masks moved instead of deep-copied when the list grows
   (it was quadratic in 300 MB copies); piece sizes counted in one pass, not one full-box count per piece;
   Dice from cached voxel lists, not full-box counts per pair; first-level regions in one label pass.
3. Threading: distance_transform lines run in parallel (exact); candidate overlap sums per slice (deterministic).
4. Candidate masks computed in the core's box padded by MDIL + 8 MFALL (values dropped < 1e-13).
5. Site layer: the single-job site rankings (fs_sync ranks, pooling climbs) are capped by a core budget
   nodepar x nparts x nthr (256 cores for the NOT keys) instead of nodepar jobs. Multi-part steps keep nodepar.
   Found because the per-site layer ran the NOT root/r2 sites 4 x 8 cores at a time for hours (the median-32-cores tail).
6. Fix found while testing: srate(:) was sized S_SBATCH; the new step index needs S_SDISC (segfault otherwise).

## End-to-end: NOT node r2 (37,679 rows) as a maxdepth-0 tree, NOT keys (site_layer/site_states, nparts 8, nthr 8, nodepar 4)
Same fixture, seed 1, CONTRAST=maps; run dirs dist_sdisc/fix_r2_{old,new5,new6}; cores from sacct over each run's worker jobs.

| arm | build | to node leaves | whole run | worker core-h | mean cores | peak | idle >2 min | master MaxRSS |
|-----|-------|----------------|-----------|---------------|------------|------|-------------|---------------|
| old  | build_guard2 (production) | 5.54 h | 6.49 h | 353+ | 58 | 544 | 1.00 h | 54 GB |
| new5 | discovery in a worker + faster discovery | - | 5.97 h | 313+ | 66 | 544 | 0.20 h | - |
| new6 | new5 + site-layer core budget | 2.07 h | 3.08 h | 358 | 116 | 592 | 0.20 h | 9.3 GB |

Correctness: sites.txt byte-identical in all three (85 sites); old vs new6 leaf_r_s{1,2,3}_rows.txt byte-identical
(18068 / 11366 / 7575), same strongest-pole margins (288.819 / 28.501 / 52.317). new6 site layer: 2 of 85 sites kept.
Site-layer verdict (1_hetero_tree/hetero_tree_sites.txt): md5-identical in old, new5 and new6; all three wrote
STATES_DELIVERED. Whole-run wall: old 6.49 h, new5 5.97 h, new6 3.08 h.

Peak 592 cores comes from the existing node-level batch (all 85 site ranks at 16 cores submitted at once, unchanged code)
overlapping the site-layer ranks for ~10 min; it fell under 600 by itself.

## Correctness on RNase P (closed_r11_routeA fixture, maxdepth 0, site_layer/site_states on, CONTRAST=maps)
production build_guard2 vs build5: sites.txt and all four leaf_r_s*_rows.txt byte-identical (md5), leaves
38259/33682/32533/32169, same strongest-pole margins (8.381 / 3.201), 3 of 5 sites kept. Wall 1h48 vs 1h51
(discovery is small at this box, so no speed gain expected there).

## Recommendation for the next full tree run
- Build: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/dist_sdisc/build6/install (21fac34db + the worktree diff,
  x6342 build, fast gate 15/15). Keys unchanged (nodepar=4 now gives the site layer 256 cores; leave it).
- Master job: 2 cores, 24 GB is enough (r2-size master peaked at 9.3 GB; site discovery now asks 64 GB in its own job).
- Expected on the NOT tree: the 2.25 h of zero-worker gaps go, sibling nodes no longer wait for a discovery, and the
  site-layer tail runs ~7x wider. On the r2 fixture this halved the wall (6.5 h -> 3.1 h) at the same core-h, so the
  12-13 h NOT tree should land near 5-6 h; not yet measured on a full tree.
- Not done (item 3, batching short steps): the per-site jobs are 8 min, not the bottleneck any more; worth revisiting
  only if a full run still shows a long tail.

## Afternoon 2026-11-01: capping the peak, then per-step profiling (builds 7-9, uncommitted)

The NOT tree 54765571 (original code, 133,754 root particles) is at 24.9 h, still running: workers ran 12.5 h,
the master ran alone (old in-process site discovery) for 12.4 h. Site ranks are ~1,000 of its ~1,200 worker core-h.

Changes on top of build6:
7. build7: the node-level site batch (all sites ranked at once) shares the core budget nodepar x nparts x nthr with the
   site layer: S_SBATCH submits what fits and tops up as jobs finish (sb_fill / busy_update; the budget is refreshed
   from each submitted dir's status every poll).
8. build8: site ranks ask 1 thread per 5,000 particles (was 2,500; reruns: 16 threads no faster than 8 at 37k) and the
   site layer reuses the node batch's ranking whenever every particle is fitted (junk_screen=no no longer required).
9. build9: site_em, exact up to rounding: the all-particle maps are summed once per EM (they were re-summed every
   iteration) and class 1 is all - class 2 (maps are linear in the weights): 2 scatter passes per iteration instead
   of 6; the N* grid's nested top-N maps are grown, not rebuilt.
   (build10 tried compact per-thread accumulators: identical output, no speed gain; reverted.)

NOT r2 fixture (same keys, seed, CONTRAST=maps):

| arm | wall | peak cores | worker core-h | labels / merged / site verdicts / unplaced |
|-----|------|------------|---------------|--------------------------------------------|
| old    | 6.49 h | 544 | 353+ | reference |
| build6 | 3.08 h | 592 | 358 | md5-identical |
| build7 | 3.30 h | 272 | 352 | md5-identical to build6 |
| build8 | 2.53 h | 264 | 165 | md5-identical to build6 (leaves: same 3, same populations) |
| build9 | 2.26 h | 376 | 143 | md5-identical to build6 (leaves: same 3, same populations; STATES_DELIVERED) |
| build12 (pushed 9005a8c77) | 1.92 h | 256 | 121 | md5-identical to build6 (same 3 leaves; STATES_DELIVERED) |

Build9's peak: 9 site-layer class-map reconstructions (pre_rec, 40 cores each) at once. Those distributed steps were never
counted against nodepar (fs_nrun skips them, old code too: build8 had 6); faster ranks let more overlap.
11. build11: both class-map submissions (pre_rec, refine_rec) go out only while the core budget has room, or when
    nothing runs.
12. build12 (user decision): junk_screen drops the report-only identity test (20 shuffles), rotational-smear fit and
    per-member evidence gap; site_junk.txt loses their two columns. Reruns vs build9: seeds, poles, posteriors, tops and
    every EM iteration count identical, junk flags (columns 1-2) md5-identical; junk screen ~210 s -> 4 s, site rank
    748 -> 561 s and 635 -> 480 s. Fast gate 15/15.
PUSHED (user ask) to fork/hetero-tree 21fac34db..9005a8c77 as three commits (788b26b46 distribution + core budget,
6afd4dc14 site_em, 9005a8c77 junk cut); build12 == 9005a8c77 file for file.
NOT r2 fixture on build12 (fix_r2_new12, = pushed 9005a8c77): 1.92 h wall, peak 256 cores (the budget, at 0.7 h),
120.7 worker core-h; labels / merged / site verdicts / unplaced md5-identical to build6; same 3 leaves 18068/11366/7575;
STATES_DELIVERED. Versus production (old): 6.49 h -> 1.92 h, peak 544 -> 256, core-h 353 -> 121.

RNase P closed_r11 replay (hard gate) on build9: hetero_tree_labels / sites / merged md5-identical to build5 (which is
md5-identical to production), same 4 leaves 38259/33682/32533/32169: closed state retained. Wall 1h07 (build5 1h51).

Build6's node chain lost 12 min to a queue wait at fflex (cluster load, not code; absent in build7/8).

site_rank profile (one rank re-run alone, timestamped log lines, 37,679 particles, 8 threads):
- ~50 % EM iterations (0.51 s/iteration at build7, 0.34 s at build9), ~27 % one junk_screen, the rest load,
  contrast fit and the N* grid. 16 threads: EM 0.26 s/it but the junk screen slows down, so 8 threads is the best.
- build7 -> build9 on the two heaviest real sites (nsite_r_s2 1,060 iterations; the seed site 782):
  956 -> 748 s and 770 -> 635 s. Seeds, poles, pole rule, contrast, starts and every EM's iteration count identical;
  posteriors within 1e-5 (the printed digit). One noise-level non-pole endpoint (z3lo, holdout correlation ~0.02)
  picked a different N* on the seed site; non-pole tops are never read by the tree.
- The junk screen's cost (~200 s per rank where it runs: 41 of the fixture's 181 ranks, including the critical-path
  pre and refine ranks) is the identity test (20 shuffles) + rotational-smear fit, which are report-only since
  91d298030 / d09352995; the tree reads only the zero-contrast junk flag. Skipping them in the tree would change no
  tree result; left for the user to decide.

## Follow-ups
- pad2() site names past 99 (see above).
- The node-level "all sites ranked at once" batch is uncapped (85-114 jobs x 16 cores); fine on csbdevel/norm today,
  but it sets the peak.

## Measurements: site discovery on the NOT tree's own recB0 maps (box 420, 8-thread x6342/x6740P job)
Same inputs as the master used (node recB0 half-maps + b.simple states); compared with the master's sites/.

| node | master in-process | build1 (worker, same code) | build3 (exact speedups) | build4/5 (+ box-local candidates) | peak RSS |
|------|------------------|----------------------------|-------------------------|-----------------------------------|----------|
| r23 (48 sites)  | 884 s  | 600 s  | 493 s  | 144 s | 18 GB |
| r2 (114 sites)  | ~3370 s (56 min idle gap) | 2462 s | 1358 s | 346 s | 40 GB |

Correctness: build1-3 sites.txt and every site mask bit-identical to the master's (r23 48/48, r2 99 masks
compared, 15 not on disk, see pad2 note). build4: same sites, same rows, max |mask diff| 3.3e-15, 0 voxels change side of 0.5.

Bug found (not fixed, pre-existing): pad2() is write(I2.2), so every site >= 100 is written to
node/sites/site_**.mrc: in the NOT run node r2's sites 100-114 were all ranked on site 114's mask. Fix = wider
padding where site masks are named (changes results for nodes with >= 100 sites); left as a follow-up.
