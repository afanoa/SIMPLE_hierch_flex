---
name: split-test-accepts-quality-split-on-null-2026-10-27
description: "KEY RESULT: on the qtoy2-matched NULL (one conformer) the node procedure reproducibly splits by particle QUALITY (SNR) and the published state_split ACCEPTS it (K*=2); its shuffle null is view x defocus stratified, not quality; candidate fixes SIMPLE_SPLIT_CONTRAST=yes / quality strata"
metadata:
  type: project
---

2026-10-27, start-point session on routes/toy55/simqtoynull (60k, conf01 only, qtoy2 quality model):
- 8 random-K2 greedy restarts (A3+B3): 7/8 converge to one split, state SNR 1.25 vs 0.82, NMI with SNR quintiles 0.19-0.21,
  junk ~0.1 %. Confidence 0.86-0.90 (higher than real r1211 0.81 / OLE 0.75); halves gate 0.96-0.98 (same as real).
- Published state_split on it: pair 1.00 vs thr 0.230, cross-half EM re-test 0.66 vs 0.096 -> K* = 2 (30,610 / 29,390).
- On qtoy2 itself (structure + same quality spread) refinement finds structure (SNR NMI 0.000), even from quality-split starts:
  quality is the fixed point only when no structural signal is left -> deep / small nodes are where it bites (fits RNase P r3
  and the m33 root site 6 being quality axes).
- Truth-free readouts (restart confidence, cross-restart reproducibility, halves gate) cannot tell it from real heterogeneity.
Runs: start_point_v1/runs/qnull_R*_g_s1; scores start_point_v1/scores/score_t6.json.
**How to apply:** the split test needs a quality-aware null before any stop/merge rule is trusted; candidates (owned by the
adaptive split/merge session): SIMPLE_SPLIT_CONTRAST=maps|yes (earlier homogeneous control: =yes merged it, some purity cost),
or a per-particle quality stratum in the shuffle null. Related: [[null-control-fixtures-2026-10-27]],
[[tree-design-feature-freeze-2026-10-25]], [[full692-tree-vs-deposited-benchmark-2026-10-22]].
- FIX FOUND (adaptive session, 49 state_split runs, ole_12707/split_merge_replay/qtest/): SIMPLE_SPLIT_CONTRAST=yes (per-particle
  amplitude in the split test) rejects the pure-SNR null splits (R1, R8: K*=1) and keeps all 4 real splits (qtoy2 root, qtoy_H,
  RNase P r12111 closed-enriched, m33 closed pair) on m8c90t and the integration source. contrast=maps does nothing; a corr-tertile
  quality stratum adds nothing on top. Residual null R4 (48,348 / 11,652) survives contrast=yes but is NOT quality (snrq AUC 0.46,
  bfac 0.50): some other nuisance (view/defocus/pose error?). Next: full trees with contrast=yes on nulls, qtoy2, OLE + equiv merge-up.
- FULL-TREE A/B on intB (integration + sequential null, seed 1, site_layer=no; adaptive session, ole_12707/split_merge_replay/trees/):
  OLE null + SIMPLE_SPLIT_CONTRAST=yes -> ONE leaf (root re-test 0.05 vs thr 0.106); unset arm splits the root 42,867 / 37,732
  (T 0.81 vs 0.099). qtoy null + contrast=yes still splits the root, weakly (T 0.09 vs 0.064; unset 0.27 vs 0.066).
- qtoynull_base_s1 (d3513ff70, unset): 5 leaves; today's merge-up 4 units; equivalence merge-up keeps all 5 (rule only removes
  merges, the false structure is in the splits: root z 15.6, r1 5.8, r2 3.4).
- RESIDUAL null split under contrast=yes (qtoy null root 31,841/28,159, z ~5) is FEATURELESS: AUC 0.46-0.53 for view,
  defocus, true pose error, shift, corr, snrq, bfac. Suspected leak (state_mergeup%test ~l.684-716 on 8c1b66f31): the
  re-test's shuffle null scores shuffled stage-B labels directly and never re-runs the cross-half EM, while the real labels
  are EM-optimized on shells 2..12/17 and re-tested 3 shells higher inside a small focus mask (mask couples ~box/d shells) ->
  anti-conservative (the 8/pi selection effect). Decisive test (adaptive session): EM-inclusive shuffle null on the null root
  and on positives (qtoy2 root, r12111, m33 closed pair).
- 16:00 adaptive session: EM-inside-the-null is degenerate (shuffled starts drain to one group, positives too); the leak is RULED
  OUT (cross_em is truly cross-half; high-pass and per-shell orthogonalisation leave T 0.09). The qtoy null splits are real but
  tiny differences (misaligned-particle share 3.7 % vs 6.3 %; R4 nothing): with ~30k/group the test answers "any difference".
  Candidate fix: relative effect rho = per-particle difference SNR / consensus SNR in the re-test band, floor 0.03 (env
  SIMPLE_SPLIT_RHOMIN). Nulls 0.012-0.023; positives m33 closed pair 0.046 (thin, 1.5x), r12111 0.116, qtoy2 root 0.555,
  qtoy truth splits 0.95-1.47, OLE r12111 1.33 / r12211 0.45. Full trees contrast=yes + rho 0.03 on qtoynull/olenull/qtoy2/OLE
  running (trees/<ds>_fix_s1); rho on all 23 m33 distinct-sibling nodes to price the RNase P cost.
- rho on all 23 RNase P m33 distinct-sibling nodes (contrast=yes build): 17 pass the re-test with rho 0.039-0.187 (closed pair
  0.046) -> a 0.03 floor costs NO RNase P split; across 3 nulls / 32 real: nulls <= 0.023, real >= 0.039 (~1.3x margin each side;
  seed stability at the margin untested). The other 6: 3 K*=1 at the board (junk screen), 3 (r1211112, r2111112, r2121212) now
  FAIL the re-test (z < 1.67) on this build though they split on m8c33/36 -> check whether they are real (deposited matches).
- The 3 m33 splits that fail the re-test on the test build (r1211112 7.3k/4.8k, r2111112 5.0k/3.4k, r2121212 3.2k/3.2k) are REAL
  and UNIQUE (amp_mstep_v1/motion_cover.json, raw unit maps 15 A): sibling difference half-reproducibility 0.65-0.70; no other
  pair of the 81 leaves carries the same difference (best 0.16-0.28, median 0.04); siblings are each other's nearest state
  (dev corr 0.71-0.78 vs cousins 0.56-0.68). Maps delivered: ~/Downloads/rnp_m33_<leaf>_unit<NN>.mrc (emp2, 10 A). Pending:
  contrast unset vs yes on these nodes to find what loses them.
- DECISIVE (17:xx): SIMPLE_SPLIT_CONTRAST=yes is what kills the 3 real unique m33 splits (same saved labels, deterministic):
  m8c36 / d3513ff70 / test build unset all K*=2 (unset rho 0.111 / 0.248 / 0.137); contrast=yes K*=1 (rho 0.000 / 0.089 / 0.042).
  Likely the per-particle amplitude fit absorbs real occupancy-like differences inside the focus mask. The qtoy null quality
  split without contrast has rho ~0.1, so a floor alone doesn't stop it. Next (adaptive session, test build only): contrast
  estimated OUTSIDE the focus region, + rho floor; battery of 31 cases then full trees.
- rho with contrast UNSET (34 cases): nulls split with LARGE effects (OLE null root rho 0.190, qtoy null root 0.104, R4 0.049)
  while RNase P real splits are 0.039 (closed pair, the lowest) to 0.255 -> no floor separates them without contrast handling.
  So: unset fails the null, contrast=yes kills real RNase P splits; candidates left = ampg (whole-molecule amplitude, fixed per
  particle; SIMPLE_SPLIT_AMPGLOBAL, ~/SIMPLE_smr5g) and maps-only; decision table pending.
- DECISION TABLE (K*/rho; unset | maps | ampg | yes): OLE null root 2/.190 | 1/.011 | 1/.011 | 1/.012; qtoy null root 2/.104 |
  2/.022 | 2/.021 | 2/.023 (floor refuses all contrast-aware arms); closed pair kept in all (.039-.046); r12111 kept except ampg;
  unique m33 splits: r1211112 kept unset/maps, r2111112 kept unset/maps/ampg, r2121212 kept ONLY unset (maps z 3.50 -> 0.93).
  BEST = SIMPLE_SPLIT_CONTRAST=maps + rho floor 0.03: rejects both contrast-carrying null roots, keeps closed pair, r12111,
  2/3 unique splits, all qtoy2/OLE positives; loses r2121212. ampg worse than maps; yes loses all 3 unique. Full trees
  trees/<ds>_fixm_s1 (qtoynull, olenull, r1211, qtoy2, OLE). r2121212 contrast-dependence check requested.
- r2121212 is NOT a contrast split (contrast AUC 0.497, corr 0.52, defocus/views ~0.5; difference reproduces 0.78 with and
  0.78 without per-particle contrast, maps correlate 0.999). The maps arm loses it through POWER: near-identical EM labels, but
  the protected re-test band (shells 16..23 = 14.6-10.2 A at box 160) holds little of a difference that lives at <= 15 A, so
  T 0.11 -> 0.05 (z 3.5 -> 0.9) on a small label change. Seed-stability battery (SIMPLE_SPLIT_SEED, 90 runs) running; if the
  decisions are coin tosses the fix is power (re-test band / replicated EM), not the contrast mode.
- Tree level (20:30): seeds 10/10 = power limit, CONFIX hypothesis wrong (E-step already a*CTF*y). OLE null tree: unset 6 leaves
  (largest 35 %), contrast=yes+rho0.03 and maps+rho0.03 both 1 leaf 100 %. maps+rho0.03 agreed; r2121212 loss accepted.
- 21:00 BOTH nulls pass at tree level on maps+rho0.03: olenull_fixm 1 leaf 100 %, qtoynull_fixm 1 leaf 100 % (root rho 0.022 < 0.03).
  Pending: r1211/qtoy2/OLE fixm trees (closed gate, real splits). Next: median-over-EM-restarts switch, battery on saved nodes first.
- 22:00 median over 4 EM restarts REJECTED: r2121212 not rescued, loses r1211112/r12111, +40-100 % split-step time. Standalone
  state_split uses ntests 1 (z_crit 1.67); in-tree z_crit ~3.1-3.4, so r2121212 fails in-tree even unset. Only (3) moving the
  protection boundary down could recover it; needs user go-ahead.
- 23:00 qtoy2 full tree (maxdepth 7) on maps+rho0.03: 18 leaves conformer NMI 0.606 (9 units 0.663) vs 0.51 base (depth 1),
  0.44 older depth-matched; purer best leaves for 6 of 8 conformers; both floor refusals (rho 0.022/0.024) had NO conformer
  signal (NMI <= 0.007). r1211 closed F1 .44 vs .45 base. Remaining: OLE fixm landscape.
- 2026-10-28 02:00: at 3+3 (recipe4 = fixm + 3+3 + crops + embed fix) the qtoy null FAILS: root K*=2 (45,494/14,506, largest
  76 % < 80 % bar), rho = 24.865e-6/564.706e-6 = 0.044 > 0.03 floor. OLE null passes (rho 0.011, 1 leaf). The fixed floor is not
  null-safe at 3+3; quality-aware pseudo-split null (SIMPLE_SPLIT_QNULL, adaptive session) is the candidate fix.
- 2026-10-28 02:30 TRAP: split Bonferroni count = 2*N_root/minpop - 1 (692k -> 460 tests, z_crit 3.71). With maps contrast the m33
  closed pair r1111221 re-tests at z 3.66 -> FAILS in a 692k tree (the r1211 fixture's count is only 18). Quality null refuses both
  null roots (qtoy null split = alignment-quality ordering, snrq AUC 0.50). Overnight 2x2 (contrast unset/maps x qnull off/on) +
  node-local count alternative running; user decides.
