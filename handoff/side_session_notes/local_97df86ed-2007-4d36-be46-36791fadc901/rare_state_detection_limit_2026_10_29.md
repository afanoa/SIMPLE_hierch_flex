# Rare discrete state: hetero_tree detection limit (session tree-detection-limit, 2026-10-29)

## Verdict
- **Detection statistic:** n / sqrt(N) (rare count over the square root of the node size), not the count or the fraction alone.
  For an ISOLATED discrete state with the closed-arm effect (conf18 vs conf01, rel. difference 0.60 at 40-20 A), at qtoy SNR:
  - found at n / sqrt(N) >= 15.5 (20 % of 6k; 10 % of 60k);
  - lost at the root Bonferroni bar at 11-12 (10 % of 12k; 5 % of 60k is at the bar, z 3.42, with only weak enrichment);
  - lost at discovery at <= 7.3 (every fraction <= 4 % at 30-60k, 3 % at 60k).
  - Extrapolated to 692k: ~1.4-1.9 % (n ~ 10-13k). The real closed state (1.26 %, 8,720) is just below that as an isolated
    state, consistent with the base tree losing it and route A / hierarchy enrichment being what saves it.
  - Smaller effects need far more: the half-closed arm (conf22, 0.43) is lost at 15.5 and the adjacent arm (conf02, 0.21) at
    24.5, both at discovery.
- **Count vs fraction:** at a fixed count of 1200, raising the fraction (2 -> 4 -> 10 -> 20 %) goes lost -> lost -> bar -> found.
  At a fixed fraction of 10 %, raising N 5x (12k -> 60k) goes bar -> found. Both matter, combined as n / sqrt(N).
- **Losing step:** discovery below ~10 (flex/fflex latent does not rank the rare particles; no site EM start or pole reaches
  them; route A has no candidate; the tree falls back to "the two groups of its fflex states", which carry the rare state only
  once flex ranks it). At ~11-12 stage B concentrates it (2.3x base) but the root re-test misses z_crit 3.71. With that bar
  removed, the locality guard is next: a diluted rare split reads L/L_ref 0.40-0.49 (merged), against quality splits at 0.20-0.34
  and found splits at 1.26. Above ~15 it passes all three.
  - minpop, choose routing and merge-up never decided a loss. The route A minpop rule is not the bottleneck either: no route A
    pole formed at any point, and the found cases were found WITHOUT route A, via the fflex fallback and enrichment.
- **No clear input-free lever.** Hierarchical Bonferroni (SIMPLE_SPLIT_HIER=yes, node at N_root / N_node tests;
  hier_ntests_on_r17d.diff, md5 482267a5, default off, build detlim/hier/build_hier, fast gate 15/15):
  - L10_12k now passes the bar, but the guard merges it (0.397); L20_6k unchanged.
  - All four nulls stay 1 leaf: qtoynull (z 5.45, guard 0.205), olenull (4.74, 0.224), N00_12k (2.17, 0.246), plus N000 at
    production keys.
  - site_dirs = 1..10 (site_rank replay on L10_12k's two arm sites): no reproducible endpoint (hold-out 0.04-0.19 vs SE 0.09).
  - The only knob with a gap is the guard threshold: quality <= 0.34 vs diluted rare >= 0.40, a 0.06 margin on few nulls.
    Not proposed. No default change.
- Figure: rare_state_detection_curve_2026_10_29.png (left: every run by n / sqrt(N) and effect size; right: conf18 in fraction x
  count with the n / sqrt(N) = 11 and 15.5 lines).

Work dir (batch2): /mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim/ (scripts/, pool/, fix/, runs/<arm>/<fixture>, traces/<arm>/, logs/).
Build: recipe17d (build_r17d/install, fast gate 0). Env: SIMPLE_SPLIT_CONTRAST=maps only (RHOMIN/ROUTE/LOCALITY unset). Launcher sst_gate4.

## Fixture (reusable)
- pool/ = one 63,600-row qtoy-model simulation (scripts/sim_pool.sbatch, job 54663525): conf01 x 59,976 majority, 1,200 each of
  conf18 (closed arm), conf22 (half-closed) and conf02 (adjacent arm position), qtoynull's junk (12 blob + 12 broken). Real RNase P poses
  and CTF (make_oris_pool.py, seed 17), noise-free simulate_particles per class, rare_merge.py = qtoy_merge's model (B ~ U(0,200),
  q lognormal 0.35, SNR0 0.0012) plus a cidx column; refine3D (prob, lp 8, 4 its) from 10 deg / 3 px perturbed truth, ONCE on the pool.
- fix/<TAG>/ (make_fixtures.py) = sst_gate4 FIX layout (input.simple = the pool project, fixture_selection.txt, fixture_rows.txt,
  truth.txt, manifest.txt). Each fixture is 60,000 rows: junk + the first (59,976 - n) majority + the first n of one rare conformer,
  so smaller rare sets nest inside larger ones. L = conf18, M = conf22, S = conf02; 200/100/050/025 = 2 / 1 / 0.5 / 0.25 %
  (n = 1200/600/300/150); N000 = no rare state (the null at the same keys).
- Effect size vs conf01 (scripts/effect.py, rel. difference norm in a 232 A sphere): conf18 0.60 (40-20 A) / 0.61 (40-10 A);
  conf22 0.43 / 0.49; conf02 0.21 / 0.30. conf02 carries ~1/8 of conf18's low-resolution signal energy.
- Caveat: the pool is refined once, so every fixture's input poses come from a refinement whose reference holds ~1.9 % of each of the
  three rare conformers. Pose quality is the same for all fixtures; detection is decided by the tree's own per-node refinements.
- Launch: scripts/launch_battery.sh (MD, MINPOP, SEED, ARM). Score: scripts/trace_runs.sh ARM TAG... -> traces/<ARM>/<TAG>.txt
  (rare_trace.py: per node flex AUC / tails, every site's seeds + pole rule, route A extra seeds, seed / stage-A / stage-B / split
  groups vs truth, re-test and guard lines; leaves with purity / recall / F1 and the map readout; merge-up units).

## Operating point
minpop = 260 = production's minpop fraction (3000 / 692,392 = 0.43 %), which also reproduces production's Bonferroni count at 60k
(2 N / minpop - 1 = 460, z_crit 3.71). qtoy's gate minpop 900 (1.5 %) would cut every rare set below 1.5 % at route A alone.

## Code facts that bound the limit (from reading recipe17d source)
- route A extra seed (state_seeds): pole equivalent-pure-members >= minpop, Spearman with the seeding site < 0.5, extra states < 1/2 node,
  driving seeds keep >= minpop.
- site_rank: starts = latent tails of min(15000, n/6) particles; N* grid SITE_GRID starts at 500 (no pole smaller than 500 members);
  pole rule: share < 0.5 and halves accept z >= 3 (row and column).
- discover_sites: sites from the variance between the K = 8 fflex state maps (count-weighted); a rare state shows only if an fflex
  state concentrates it.
- ntests = 2 N_root / minpop - 1; decide(): children only when every split group >= minpop (else the groups become leaves).

## Results: waves 1-2 (maxdepth 1, seed 1; runs/<arm>/<tag>, traces/<arm>/<tag>.txt)
Every run below ends as ONE leaf (100 %); the loss is root discovery in all of them.
| tag | N | rare (n, %) | flex best AUC | sites / poles reaching rare | route A | stage B rare share s1 / s2 | re-test (z_crit 3.71) | guard | leaf F1 | core-h |
|---|---|---|---|---|---|---|---|---|---|---|
| N000 | 60k | none | - | - | - | - | - | - | 1 leaf (null passes) | 11.6 |
| L025 | 60k | conf18 150 (0.25) | 0.61 | none | none | 0.2 / 0.3 % | PASS 9.70 | 2 -> 1 | 0.005 | 14.4 |
| L050 | 60k | conf18 300 (0.5) | 0.41 | none | none | 0.5 / 0.5 % | FAIL 2.80 | - | 0.010 | 12.5 |
| L100 | 60k | conf18 600 (1) | 0.41 | none | none | 0.8 / 1.3 % | PASS 7.06 | 2 -> 1 | 0.020 | 11.4 |
| L200 | 60k | conf18 1200 (2) | 0.40 | none | none | 1.5 / 2.6 % | PASS 4.38 | 2 -> 1 | 0.039 | 11.7 |
| S200 | 60k | conf02 1200 (2) | 0.48 | none | none | 1.8 / 2.2 % | PASS 5.24 | 2 -> 1 | 0.039 | 12.0 |
| H200 | 30k (minpop 130) | conf18 600 (2) | 0.39 | none | none | 1.8 / 2.2 % | PASS 6.58 | 2 -> 1 | 0.039 | 7.8 |
| H400 | 30k (minpop 130) | conf18 1200 (4) | 0.42 | none | none | 3.4 / 4.6 % | FAIL 3.56 | - | 0.077 | 5.9 |
- Path in every run: flex/fflex latent does not rank the rare particles (no eigenvalue spike vs the null; no fflex state above
  ~1.4x the base purity); every site's starts score AUC ~0.5 for them and no pole holds them; route A has no candidate pole; the
  node falls back to "no site gives a split, the two groups of its fflex states seed the refinement"; stage B refines a quality-like
  2-split; the re-test refuses it or the locality guard merges it (correct: the two maps differ diffusely). minpop, z_crit, choose
  routing and merge-up never see the rare state.
- Wall 19-42 min per run, ~12 core-h at 60k, ~7 at 30k (K* = 1 roots, no children).
- Reading: an ISOLATED discrete state (no correlated axis) is below root discovery at <= 4 % at qtoy SNR, at 30k and 60k. qtoy2's
  6 % closed state is found because it sits at the end of the arm-swing axis most particles share (root z3 AUC 0.96, fflex state 8
  purity 0.87); the real 692k closed state (1.26 %) is found after enrichment down the hierarchy. Spiked-covariance estimate: the
  rare state's spike f x d'^2 (d' ~ 0.3-0.7 per particle at this SNR) sits far below the noise bulk edge ~ sqrt(p / N).

## Results: wave 3 (fixed count 1200, smaller N; minpop 0.43 % of N -> 460 tests; 29.5 core-h)
| tag | N | rare | n/sqrt(N) | flex best AUC | stage B best state | re-test | outcome / losing step |
|---|---|---|---|---|---|---|---|
| L05_12k | 12k | conf18 600 (5 %) | 5.5 | 0.67 | 5.5 / 4.3 % | FAIL 1.86 | 1 leaf; discovery |
| L10_12k | 12k | conf18 1200 (10 %) | 11 | 0.30 (z4) | purity 0.23, recall 0.91 | FAIL 3.29 < 3.71 | 1 leaf; BONFERRONI bar |
| L20_6k | 6k | conf18 1200 (20 %) | 15.5 | 0.12 (z1) | purity 0.57, recall 0.95 | PASS 5.97, choose -> EM groups, guard keeps | FOUND: leaf r2_g1_2 1449, P 0.71 R 0.86 F1 0.78, map rare 0.74 / maj 0.47 |
| S10_12k | 12k | conf02 1200 (10 %) | 11 | 0.53 | 10 / 10 % | FAIL 0.27 | 1 leaf; discovery |
| S20_6k | 6k | conf02 1200 (20 %) | 15.5 | 0.48 | 21 / 19 % | FAIL 0.34 | 1 leaf; discovery |
| N00_12k | 12k | none | - | - | - | - | 1 leaf (null passes) |
- L20_6k path: no site pole kept (the arm site exists: its mask holds 22 % of the rare-difference energy in 0.7 % of the box), the
  fflex-group fallback carries 84 % of rare (purity 0.29), stage A/B purify (0.48 -> 0.57). Child r2 (71 % rare) re-test z -1.51:
  it cannot purify further; 165 rare stay in the majority leaf.
- L10_12k: sites 3-5 sit on the arm (17-41 % of the rare energy) but every site-EM start (z1 / z3 tails, 2000 each) gives a
  non-reproducible endpoint (hold-out corr 0.04-0.11, halves accept z -1.9 .. 1.8 vs 3), so no pole; the rare direction is z4 / z7,
  which site_dirs=1,3 never starts from. The fallback carries 69 % into stage B, which purifies to 0.23 / 0.91.
- Replay (replays/L10_12k_r_m*, state_split on the saved stage B, guard skipped: maps pruned): m 460 z 3.29 FAIL (= tree);
  m 132 z 3.29 vs 3.38 FAIL; m 1 z 3.57 vs 1.67 PASS -> K* 2: group 4681 rows, rare 1088 (P 0.23, R 0.91).
- Count vs fraction: at fixed count 1200, 2 -> 4 -> 10 -> 20 % goes lost -> lost -> bar -> found, so count alone does not decide.
  All points fit n / sqrt(N): lost <= 6.9, bar at 11, found at 15.5. Extrapolation: ~1.3-1.8 % at 692k for an isolated state
  of the closed-arm effect. Fixed-fraction test at 60k (pool2: 10 / 5 / 3 % conf18, 10 % conf02) launched 19:30 (job 54680482).
- Small effect (conf02, 1/8 of the low-res signal energy): flex never ranks it up to 20 % (n/sqrt N 15.5): discovery.
- TRAP: hetero_tree's end-of-run prune deletes split/s.simple and every stage-B map; replays then fall back to stageB/b.simple and the
  locality guard cannot run. Launch with SIMPLE_TREE_KEEPMAPS=yes (launch_battery.sh now exports it).

## Results: child cascade, lever, count test, medium effect (maps kept from here on)
- C_L10m1_g2 (L10_12k's replayed m=1 group 2 as a fixture, 4681 rows, 1088 rare, split_ntests=3): fflex fallback -> stage B
  purity 0.41 / recall 0.79 -> re-test PASS z 2.39 vs 2.15 -> choose EM groups 2972 / 1709 -> guard L/L_ref 0.487 -> merged, 1 leaf.
- Lever runs (build_hier, SIMPLE_SPLIT_HIER=yes, maxdepth 2; runs/hier_*):
  - L10_12k: root z 3.57 vs 1.67 PASS -> guard 0.397 -> merged, 1 leaf.
  - L20_6k: identical to production (found; child re-tests at m 1 / 4 FAIL).
  - N00_12k z 2.17 -> 0.246 merged; qtoynull z 5.45 -> 0.205; olenull z 4.74 -> 0.224: every null 1 leaf (100 %).
- site_rank replay (replays/L10_12k_site{3,4}_all, site_dirs 1..10, 20 starts per site): every endpoint rejected (no pole);
  the control site3 with 1,3 reproduces the tree.
- pool2 (conf01 59,976 + conf18 6,000 + conf02 6,000 + junk, refined once, pose error 3.77 deg; job 54680482), 60k fixtures,
  production keys (minpop 260):
| tag | rare | n/sqrt(N) | flex best AUC | stage B best state | re-test | guard | outcome |
|---|---|---|---|---|---|---|---|
| P_L03 | conf18 1800 (3 %) | 7.3 | 0.42 | P 0.049 R 0.73 | PASS 4.57 | 0.201 merged | 1 leaf; discovery (quality split) |
| P_L05 | conf18 3000 (5 %) | 12.2 | 0.40 | P 0.063 R 0.57 | FAIL 3.42 | - | 1 leaf; discovery / bar edge |
| P_L10 | conf18 6000 (10 %) | 24.5 | 0.76 | P 0.25 R 0.96 | PASS 13.35 | kept | FOUND at depth 1 |
| P_S10 | conf02 6000 (10 %) | 24.5 | 0.52 | P 0.11 R 0.57 | PASS 7.47 | 0.203 merged | 1 leaf; discovery |
  - P_L10 path: root has no site pole and no route A pole; the fflex fallback carries 88 % at 14 %; stage B gives 22,922 rows at
    25 %. Depth-1 node r2: flex AUC 0.83, fflex state purity 0.83 -> stage B 7,928 rows at 64 % -> re-test z 13.21 -> leaf r2_s1
    7,928 rows, P 0.64, R 0.84, F1 0.73, map rare 0.75 / majority 0.36 (reaches). Rare elsewhere: 683 in r2_s2 (4.6 %), 259 in
    the majority leaf r1 (0.7 %).
  - The guard really evaluated every found split, so none passed through the EM-group skip trap: P_L10 r 0.712, r2 1.308,
    L20_6k 1.261. At a similar purity the guard ratio also grows with count: P_L10's root split (60k, purity 0.25) reads 0.712
    and is kept; L10_12k's (12k, purity 0.23) reads 0.397 and is merged.
- Medium effect (conf22): M10_12k (10 %, 11) flex 0.45, PASS 4.22 -> guard 0.253 -> 1 leaf; M20_6k (20 %, 15.5) flex 0.57,
  FAIL 2.39 -> 1 leaf. Both lost at discovery.
- Where the rare particles end up when lost: always in the single leaf (100 %), so the leaf F1 equals 2f / (1 + f).

## Cost and disk
- About 330 core-h: trees ~290 (waves 1-3 117, lever battery ~62, count test ~73, medium + child ~12), simulations ~60 (two
  pools). Walls 11-42 min per run.
- detlim/ is ~147 GB, over the ~50 GB norm. Mostly kept stage maps (5-14 GB per run since KEEPMAPS) and the pools.
  - Moved aside for deletion (user only): _to_delete_20261029/pool{,2}/stk_*.mrcs (noise-free per-class stacks, 14 GB; the
    fixtures read pool*/all.mrcs).
  - Also deletable once nobody needs replays: runs/*/*/1_hetero_tree/node_*/{stageA,stageB,recB,seeds,...}/*.mrc.

## Reuse
- New fixture from the pools: `N_FIX=<N> python3 scripts/make_fixtures.py pool/truth.txt pool/src.simple fix 1 TAG=CONF:n`.
  pool holds 1,200 each of conf18 / 22 / 02; pool2 holds 6,000 each of conf18 / 02.
- Run: `MD=1 MINPOP=<0.43 % of N> bash scripts/launch_battery.sh TAG ...`; nulls `ARM=x bash scripts/launch_nulls.sh`.
- Score: `bash scripts/trace_runs.sh <arm> TAG ...` (traces/<arm>/<TAG>.txt) or `null_leaves.py`.
- Replay a decision: split_replay.sbatch (state_split at any count / env) and site_replay.sbatch (site_rank at any site_dirs),
  both table-driven; replay_groups.py and child_fixture.py turn replayed groups into scores or new fixtures.
- Figure: scripts/curve.tsv + scripts/curve_fig.py.

## Correlated-axis series (rare state beyond the end of a shared swing; 2026-10-29 evening, ~168 core-h)
- pool3 (job 54692082): the qtoy2 arm swing as the majority, at qtoy2 weights and with no holo: conf43 / 49 / 28 / 01 / 02 / 22 =
  4,862 / 8,104 / 12,968 / 16,212 / 11,346 / 6,484. Plus conf18 x 1200 beyond the conf22 end and qtoynull's 24 junk. Refined once
  (pose error 3.99 deg). Fixtures A200 / A100 / A050 (2 / 1 / 0.5 % of 60k; make_fixtures MAJ 0 = all other conformers, one seeded
  row order, nested). Production env, minpop 260, maxdepth 2, maps kept (runs/mp260_md2_s1).
- Path of the rare state (share, n/sqrt(N), share of all rare):
| run | root | depth 1 (closed side) | depth 2 (end of swing) | losing step | best leaf |
|---|---|---|---|---|---|
| A200 | 2.0 %, 4.9 | r1 33.7k (01 / 02 / 22): 2.44 %, 4.5, 0.69 | r12 10.9k (22 37 %): 6.5 %, 6.8, 0.59; stage B 4,930 rows at 12.2 %, R 0.85 | bar at depth 2 (z 2.38 < 3.71) | r12_g1_2 P 0.065 R 0.59 F1 0.12 |
| A100 | 1.0 %, 2.4 | r2 16.8k: 1.77 %, 2.3, 0.50 | r22 6.8k: 3.9 %, 3.2, 0.44; stage B 5.0 / 3.0 % | discovery at depth 2 (z 1.49) | r22_g1_2 F1 0.07 |
| A050 | 0.5 %, 1.2 | r1 33.8k: 0.63 %, 1.2, 0.71 | r12 8.4k: 1.97 %, 1.8, 0.55; stage B 1.7 / 2.2 % | discovery (z 0.32) | F1 0.04 |
  - Root: the swing splits into open (28 / 49 / 43) and closed-side (01 / 02 / 22) halves (K* 3-4, guard 1.10-1.15). Route A added
    swing poles; none held conf18. Root flex ranks conf18 only moderately (AUC 0.29-0.31 at 2 / 1 %).
  - Deepened per the coordinator's rule: the A200 r12 split replayed at 1 test passes (z 2.38 vs 1.67) and the guard keeps it
    (1.125). Its group 2 (4,930 rows: conf22 2,085, 02 946, 01 722, conf18 599 = 12.2 %, n/sqrt(N) 8.5) run as a child fixture
    at the tree's count (split_ntests 460; runs/child_a): flex AUC 0.47, stage B 11.9 / 12.4 %, re-test z 0.30 -> 1 leaf. Lost
    at discovery.
- Why the shared axis does not rescue it at 60k: enrichment along the swing raises the rare share (2.0 -> 2.44 -> 6.5 -> 12.2 %),
  but n/sqrt(N) only goes 4.9 -> 4.5 -> 6.8 -> 8.5. The companions it collects are conf18's swing neighbours: conf18 vs conf22 is
  rel 0.41 at 40-20 A (effect.py vs conf22), the medium effect that the isolated series lost even at 15.5. So the
  isolated-state rule holds inside the branch, at the effect size against the branch's own majority.
- Compare qtoy2 (closed 6 % = 3,600 in swing + holo + junk at 60k): found at depth 2 (F1 0.87). In the real 692k data, closed is
  8,720 particles. If the 44k r131 node holds most of them (not checked here), n/sqrt(N) there is ~30-40, well inside the
  found regime. At 60k the same fraction cannot get there.
- Disk: the three A runs with kept maps are 30-42 GB each (109 GB); detlim/ is now ~255 GB. Nothing deleted.

## Not done / follow-ups
- One seed per point. The n/sqrt(N) boundary (11-15.5) rests on 2 found and 2 bar points; seeds 2-3 at L10_12k, P_L05 and
  L20_6k would firm it (~40 core-h).
- No 200k run. P_L10 vs L10_12k already shows the fixed-fraction count effect (5x N); a 200k at 3-4 % (predicted n/sqrt N
  13-18) would test the extrapolation toward 692k (~200 core-h incl. its pool).
- Isolated states only. A rare state riding on a correlated axis (qtoy2's closed state at the end of the arm swing) is found
  at much lower n/sqrt N, because flex ranks the shared axis.
