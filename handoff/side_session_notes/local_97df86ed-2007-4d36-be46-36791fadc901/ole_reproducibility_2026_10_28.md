# OLE reproducibility (session: OLE repro, worktree vigorous-snyder-4c4173), 2026-10-28 overnight

Write-up: /Users/atifao/SIMPLE/.claude/worktrees/vigorous-snyder-4c4173/doc/implementation_notes/planned/hetero_tree_ole_reproducibility_2026_10_28.md
Scripts: same worktree, scripts/hetero_tree_ole_repro/. Cluster: /scratch/cluster_scratch/atifao/ole_repro/ (maps, trace, replays, gates, tools, build_xd).
Nothing committed. One native diff (env switch, default off): ole_repro/tools/xd_ksel_1se.diff; built into ole_repro/build_xd
(recipe4 source) and /scratch/cluster_scratch/atifao/recipe8 (recipe7 source), fast gate 0 for both.

## What matters for other sessions
- Runs diverge ONLY at the depth-1 nodes' seeding step (site or pole), never in stage A/B or the split test: 7 of 7 depth-1
  node pairs across 4 independent seed pairs; the root never diverges.
- Three near-ties decide the seeding by noise:
  1. fflex XD K = argmax of a flat held-out curve (K 3 vs 6: 0.0002 nats/point);
  2. first-site-in-order + pole rule >= 3: a marginal pole (margins 1.2-4.8 across seeds) seeds the node in 1 seed of 4;
  3. at K = 3 the main region splits into two equal factors (margins 20-23, labels arbitrary).
- r1 replays (split ARI / stage-B map corr, same rows -> 1.5 % other rows): base .69/.87 -> .27/.52; strongest_pole .88/.98 ->
  .43/.74; SIMPLE_XD_KSEL=1se alone .44/.59 -> .35/.58; 1se + strongest .68/.88 -> .53/.75.
- 1se gates pass (RNase P r1211 closed F1 .423/.517 -> .524/.514; qtoy2 NMI .542 -> .541; K reproducible per node), but alone it
  hurts OLE r1. strongest_pole gates (r1211 s1/s2, qtoy2) running in ole_repro/gates/*_strong_*.
- Full pair running: recipe7 + strongest_pole, seeds 1/2, /mnt/beegfs/elmlund/afan/ole_repro_20261028/recipe7_strong_s{1,2}
  (nparts 2, nodepar 4, site_layer=no). Score with particle_repro.py + make_partitions.py/rec.sbatch/maps_score.py.
- Map level (consensus-pose leaf maps vs nulls shuffling below depth 1/2): span, extremes and depth-2 structure reproduce;
  deeper leaves are interchangeable samples (cross-run best match 0.64-0.82 = own nearest sibling 0.75-0.80).
- Fixed-size refined end poles reproduce only when the axis does (0.61-0.75 same axis, ~chance across independent seeds).
- Proposed next (untested): treat sibling factors of one discovered component (1f1/1f2) as one site candidate.

## Final (2026-10-28 ~11:55)
- strongest_pole alone passes every gate run: r1211 closed F1 .62/.56 (site order .42/.52); qtoy2 holo F1 0.82 (same as site
  order at maxdepth 7); olenull 1 leaf. qtoynull fails in both arms (0.669 vs 0.758). Cost: small qtoy2 conformers c43/c49 lose
  about 0.2 F1 each.
- OLE pair (a) recipe7 + strongest s1/s2 (beegfs ole_repro_20261028): depth-2 ARI 0.77 (site order 0.35-0.49); map coverage
  0.78/0.83 (15/10 A); 78-83 % of particles in states partnered >= 0.7; extremes 12/12 at 10 A; leaves 0 at J >= 0.5;
  252-270 core-h per tree (= site-layer tree).
- Sibling-factor rules (SIMPLE_SITE_SIBLINGS both/union) fail the r1211 closed gate (seed 2: 0.38/0.42): dropped.
  SIMPLE_XD_KSEL=1se kept as a test switch only.
