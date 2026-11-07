---
name: ole-repro-seeding-near-ties-2026-10-28
description: "OLE runs diverge ONLY at depth-1 seeding: three noise-decided near-ties (fflex XD K argmax on a flat curve, a marginal pole at the pole-rule >=3 threshold taken first by site order, two equal factors of one region); strongest_pole best tested fix; SIMPLE_XD_KSEL=1se alone hurts r1"
metadata:
  node_type: memory
  type: project
  originSessionId: 2f3dad73-7d79-4c2e-b1d3-038344148ed7
  modified: 2026-10-28T10:39:59.231Z
---

2026-10-28 overnight, OLE reproducibility session (worktree vigorous-snyder-4c4173; write-up
doc/implementation_notes/planned/hetero_tree_ole_reproducibility_2026_10_28.md; cluster /scratch/cluster_scratch/atifao/ole_repro/).

- Node trace on 4 independent seed pairs: the root never diverges; 7/7 depth-1 node pairs diverge at SEEDING (site or pole);
  stage A/B and the split test never start a divergence (seeds ARI >= 0.6 -> split ARI 0.79-0.98).
- Near-tie 1: fflex's XD deconvolved-mixture K (`xd_select_k`, simple_xd_gmm.f90) is an argmax over a flat held-out curve:
  K = 3 vs 6 at OLE r1 differ by 0.0002 nats/point. K -> fflex states -> stA0/recB0 maps -> discover_sites list and order.
  flex_pca draws no RNG here; tiny upstream differences decide.
- Near-tie 2: the node seeds from the FIRST site in order that keeps a pole. A marginal pole (5f1, ~6.7k rows) has pole-rule
  margins 1.2-4.8 across seeds against the >= 3 threshold, so 1 seed in 4 seeds from it (split ARI 0.05 vs the others) even
  though a margin-20-35 pole sits one site later.
- Near-tie 3: at K = 3 discover_sites factors the main 63 A region into 1f1/1f2 with equal margins (20-23); factor labels are
  arbitrary per run, so even strongest_pole tosses a coin (0.47 between groups, 0.98 within).
- r1 replays (tree 16 rows x4 seeds, tree 17 rows x2), split ARI / stage-B map corr, same rows -> other rows:
  base 0.69/0.87 -> 0.27/0.52; strongest_pole 0.88/0.98 -> 0.43/0.74; 1se alone 0.44/0.59 -> 0.35/0.58;
  1se+strongest 0.68/0.88 -> 0.53/0.75. Cost of strongest: +1.6-2.9 core-h per extra site rank (+35-47 % per OLE r1 node).
- SIMPLE_XD_KSEL=1se (env, default off; diff ole_repro/tools/xd_ksel_1se.diff; installs ole_repro/build_xd and
  /scratch/cluster_scratch/atifao/recipe8): makes per-node K reproducible (OLE and RNase P r1211); gates pass
  (r1211 closed F1 .423/.517 -> .524/.514; qtoy2 NMI .542 -> .541).
- Map level (consensus-pose leaf maps vs nulls shuffling below depth 1/2): span, extremes and depth-2 structure reproduce;
  below depth 2 leaves are interchangeable (cross-run best match 0.64-0.82 = own nearest sibling 0.75-0.80).
- Refined end poles reproduce only when the axis does (same axis 0.61-0.75, independent seeds ~chance).

**Why:** the user's goal is a reproducible conformational landscape ([[landscape-must-be-reproducible-and-splits-data-derived]]);
this localises what breaks it on OLE. Extends [[ole-seed-split-choice-2026-10-25]].
**How to apply:**
- Judge seeding changes by r1 replays on BOTH same rows and other rows (a 1.5 % row change reshuffles near-ties), and at the
  map level, not only partition ARI.
- Next candidate (untested): treat sibling factors of one discovered component as one site candidate.
- Running at write time: OLE pair recipe7 + strongest_pole (beegfs ole_repro_20261028), strongest gates in ole_repro/gates.

**Verdict 2026-10-28 (full trees + gates):**
- strongest_pole alone is the only arm passing every gate: r1211 closed F1 0.62 / 0.56 (site order 0.42 / 0.52); qtoy2 holo lost at depth 1 but recovered at depth 3 (0.97 / 0.71); olenull 1 leaf. qtoynull fails in both arms (0.669 strongest, 0.758 site order).
- Sibling-factor rules 'both' and 'union' (SIMPLE_SITE_SIBLINGS) FAIL the r1211 closed gate on seed 2 (0.38, 0.42): dropped.
- OLE pair (a) = recipe7 + strongest_pole s1/s2 (beegfs ole_repro_20261028): particle ARI d1 0.81, d2 0.77 (site order 0.35-0.49), d3 0.37, leaves 0.14, 0 leaves at J >= 0.5. Maps: extremes 12 of 12 partnered >= 0.7 at 10 A (16/17 site order 11 of 12). Strongest vs site-order trees cover each other only 0.66-0.68: the root axis differs, each method reproduces itself.
- strongest_pole's site batches burst 128-224 CPUs per node (SITE_NTHR 32 floor in code); total core-h ~ site-layer tree.
- qtoy2 like-for-like (recipe9, maxdepth 7, seed 1): site order NMI 0.615 vs strongest 0.606; holo F1 0.82 both; c18/c28 identical; strongest loses small conformers c43 (F1 0.66 -> 0.46) and c49 (0.80 -> 0.60) via its r2 cut (site 2 margin 112 vs site 1 29.7).
- Footprint: /scratch/cluster_scratch/atifao/ole_repro holds 1.4 TB of replay intermediates (not deleted).
- qtoy2 like-for-like (recipe9, maxdepth 7, site order -> strongest): NMI 0.615 -> 0.606; holo F1 0.82 both; c18, c28, c1, c2, c22 identical; small conformers c43 0.66 -> 0.46 and c49 0.80 -> 0.60 (strongest takes r2 site 2, margin 112, over site 1, 29.7). Write-up: worktree vigorous-snyder-4c4173 doc/implementation_notes/planned/hetero_tree_ole_reproducibility_2026_10_28.md; 1.4 TB replay intermediates in /scratch/cluster_scratch/atifao/ole_repro (user's to clear).
- USER DECISION 2026-10-28 ~15:20: "Adopt strongest_pole as the default; the OLE reproducibility gain is worth the qtoy2 c43/c49
  loss." -> recipe12 (optimization session local_a012d0f3) flips strongest_pole default to yes (simple_parameters.f90:65 + UI
  :787); thread-by-node-size must cover its site batches (SITE_NTHR 32 bursts).
- Choose routing (recipe14 default) vs em, OLE r1 replays seeds 1-4 (beegfs hetero_tree_20261028/choose_r1): split ARI 0.82 (min 0.66) vs 0.87 (0.76); the three like-seeded seeds 0.98 in both; only the K=3/site-3 seed loses (0.66 vs 0.76) because choose delivers stage B itself; stage-B maps 0.95-1.00 either way. Not clearly worse at the map level.
- LANDSCAPE SCORE 2026-10-28 ~19:50 (landscape_score.py, ref6 root flex_pca space, 3 axes, depth-4 states, minpop 3000; beegfs
  ole_r14/logs/landscape_54561664.out): extrema reproduce within AND across builds (max |diff| .04-.08, mean .02-.035 span units vs
  null .065-.12/.05-.07). Occupancy W1: same-build floors ref6~ref7 .051, strong s1~s2 .046; ref~recipe4 .058; ref~strongest
  .072-.080 = null range (.062-.094) -> strongest places population differently (other depth-1 axis), same ends. Coverage
  non-discriminating. Separation modest (floors 1.3-2x below null).
