# recipe18c switch matrix: results as they land (cluster: hetero_tree_20261028/recipe18c/{runs,route,score})

## ROUTE offline on the continuum (tier 2; ct_route_score.py, genuine rows; eta2 of truth t over labels, W1 of member-mean positions)
| tree | labels | eta2 | W1 | leaves (member t, n) |
|---|---|---|---|---|
| ctpk base | tree / down / flat1 / flat2 | .185 / .163 / .163 / .162 | .101 / .103 / .103 / .103 | 2 leaves .51 / .69 |
| ctpk inh | tree / down / flat1 / flat2 | .255 / .208 / .200 / .203 | .075 / .081 / .082 / .082 | 4 leaves .41 .56 .65 .72 |
| ctun base | tree / down / flat1 / flat2 | .385 / .357 / .351 / .352 | .133 / .135 / .134 / .133 | 3 leaves .31 .34 .69 |
| ctun inh | tree / down / flat1 / flat2 | .448 / .402 / .385 / .389 | .097 / .106 / .107 / .107 | 4 leaves .25 .40 .61 .77 |
- Routing moves 9-28 % of rows (own leaf .71-.91) and lowers eta2 by .02-.06 everywhere: on a continuum it blurs the tree's cuts
  (as the continuum session predicted). Cost 3.1 / 6.5 / 4.8 / 6.1 core-h.

## POSEGATE alone on the micelle fixtures (tier 1; fp inputs, md0, recipe18c)
| fixture | root stA0 flips > 30 deg | gate | re-test | K*, sizes | guard | vs fp_nomask (fixposes key, micsup_fp build) |
|---|---|---|---|---|---|---|
| micnull | 30.2 % | fired (marker, stA0_realigned) | z 2.52 < 3.01 FAIL | 1 leaf, 100 % | - | identical (z 2.52, K* 1): null PASS |
| mic3 | 28.9 % | fired | z 3.91 PASS | 2: 30066 / 26742 | 0.999 kept | identical sizes, z, guard: the micelle-size split, no protein-state separation |
- Cost 7.2 + 8.1 core-h. The gate reproduces the fixposes runs exactly; PG alone does not fix mic3's micelle split (SUPPORT does).

## Step 1: all-off recipe18c = recipe17d (L20_6k md1, identical keys, seed 1)
- Root re-test z 5.97 vs 3.71 PASS, guard L/L_ref 1.261, K* 2 = 4551 / 1449; r1 z -0.39 and r2 z -1.51 FAIL: identical lines.
- Leaves r1_g1_2 4551 / r2_g1_2 1449 (the found rare leaf); labels identical on all 6,000 rows. 8.0 core-h, 44 min.

## Nulls, SUPPORT=order + POSEGATE (tier 1; md0; baseline floor_off/*_s1 = recipe17d floor 0)
| null | gate flips | focus mask (voxels > 0.5) | re-test | guard | leaves | null bar |
|---|---|---|---|---|---|---|
| qtoynull SUP+PG | 0.08 %, not fired | 125,432 (order support, band 13.8-9.4 A) | z 1.33 < 3.38 FAIL | - | 1, 100 % | PASS |
| qtoynull base | - | 1,691,751 | z 5.45 PASS | 0.205 merged | 1, 100 % | PASS |
| olenull SUP+PG | 0.32 %, not fired | 97,442 (band 10.4-7.0 A) | z -3.77 < 3.12 FAIL | - | 1, 100 % | PASS |
| olenull base | - | 880,073 | z 4.89 PASS | 0.224 merged | 1, 100 % | PASS |
- The order support shrinks the focus mask 9-13x and the nulls' quality split no longer passes the re-test (no guard needed).
  The mask keeps only density structured at the finest third of the FSC >= 0.5 shells: the flexible-arm risk is what the RNase P
  SUP runs must show. Cost 11.4 + 28.6 core-h.

## Rare state (tier 2; detlim fixtures md1; rare_score.py; rare = conf18, 1200 rows)
| fixture | arm | root re-test (z_crit 3.71) | guard | leaves (size/rare) | best rare leaf P / R / F1 |
|---|---|---|---|---|---|
| L20_6k | base (r17d) = r18c OFF | 5.97 PASS | 1.261 kept | 1449/1035, 4551/165 | .71 / .86 / .78 |
| L20_6k | SUP | 5.75 PASS | SKIPPED (EM groups, no spanning pair) | 1532/977, 4468/223 | .64 / .81 / .72 |
| L20_6k | INH+PG (gate 0.04 %, not fired) | 5.97 PASS | 1.261 kept | identical to base (2 inherited child seedings, both re-tests fail) | .71 / .86 / .78 |
| L20_6k | ROUTE offline on base: down / flat1 / flat2 | - | - | 1038/806 ... 592/495 | .78/.67/.72, same, .84/.41/.55 |
| L10_12k | base (r17d) | 3.29 FAIL | - | 1 leaf | - (.18 = whole node) |
| L10_12k | SUP | 3.80 PASS (barely) | SKIPPED (EM groups) | 1384/843, 1188/137, 9428/220 | .61 / .70 / .65 |
- SUP finds L10_12k's rare state (base misses it) but through an UNGUARDED root split (hier-ntests run: guard 0.397 merged it);
  slightly worse at L20_6k (.72 vs .78), again with the guard skipped. ROUTE flat2 hurts the rare leaf (.55). INH no change.
- Cost: L10 SUP 9.1, L20 SUP 7.7, L20 INH 7.9, route 1.5 core-h.

## mm3q (independent motions) ROUTE offline on the recipe17d md4 tree (disc_score.py, genuine rows)
| labels | NMI(combo) | best-leaf F1 mean / min | combos >= .7 |
|---|---|---|---|
| tree (15 leaves) | .640 | .67 / .52 | 1 of 8 |
| route-down | .620 | .64 / .46 | 1 |
| flat1 | .655 | .71 / .60 | 4 |
| flat2 | .679 | .75 / .63 | 6 |
- Flat routing helps discrete independent motions (as qtoy2); route-down slightly hurts. 16.8 core-h, 41 min.

## OLE ROUTE offline on recipe17b r17b_s1 (md4; no truth)
- Own leaf after route-down 0.494, flat1 0.455, flat2 0.549; first disagreement at depth 0 for 20.8 % of rows (root K-way 5),
  depth 1 12.6 %, 2 11.7 %, 3 5.4 %. 54.5 % of rows have a min margin < 1 nat (median 1.32 agree / 0.61 disagree).
- Routing reshuffles about half of OLE's particles on sub-nat margins: no stable re-assignment to deliver. 28.8 core-h, 54 min.

## r131 (tier 1; RNase P closed subtree, md1, binary everywhere, split_ntests 460)
- Trees (cheap readout): INH+PG 4 leaves 11651/10708/11194/10781 (gate 0.33 %, not fired; 2 inherited seedings at r1, r2; re-tests
  r 19.86, r1 6.75, r2 6.15; all guarded). SUP 4 leaves 12179/10427/10961/10767 (re-tests r 12.80, r1 6.53, r2 12.27; guarded).
  Base r131_off 4 leaves. Costs 21.9 / 21.2 core-h. Map scorecard: score/r131/scorecard.txt (job 54709344).
- ROUTE offline on r131_off (route_ref.py, closed reference rows; maps not re-made): best closed-leaf F1 tree .123 / down .122 /
  flat1 .122 / flat2 .125; union .135-.137. Own leaf .69 (down), .73 (flat2); 45 % ... median margin 1.4-1.6. Neutral. 5.8 core-h.
- r131 map scorecard (leaf_scorecard.py; corner = closed-arm >= 0.95, open-arm <= 0.6):
| arm | corner leaves | best closed-leaf F1 (leaf, rows, ref, arm) | leaf arms (closed/open) |
|---|---|---|---|
| base r131_off | 1: r2_s1 0.97/0.54 | .123 (r2_s1, 10670, 845) | .97/.54 .81/.61 .81/.50 .90/.63 |
| INH+PG | 0 | .147 (r2_s1, 10781, 1017, 0.92/0.52) | .92/.52 .85/.60 .85/.64 .85/.52 |
| SUP | 0 | .125 (r2_s2, 10767, 864, 0.38/0.55) | .38/.55 .36/.56 .34/.55 .36/.54 |
- SUP FAILS the closed gate: every leaf map reads closed-arm 0.34-0.38 (the arm is cut away: the flexible-arm risk, confirmed).
- INH+PG FAILS the corner bar (0.92 < 0.95) though its closed leaf holds more reference rows (F1 .147 vs .123); one seed.

## closed_r111 (tier 1; route A r111 fixture, md1, kway 2; leaf_scorecard)
| arm | leaves | corner leaves | best closed-leaf F1 | notes |
|---|---|---|---|---|
| OFF (r18c) | 4 | 2: r1_s1 .97/.60, r1_s2 1.04/.59 | .191 (r1_s1 10688, 1388 ref) | root z 5.74 (m 18), r1 8.41, r2 7.50 |
| INH+PG | 4 | 2: identical r1 branch | .191 | gate 0.18 %; inherited 1 (r2: 10548/9720 vs 10222/10046) |
| SUP | 6 | 0 | .146 (r3_s2, arm 0.40/0.52) | root K* 3; all leaf maps closed arm 0.33-0.40: ARM CUT |
| recipe16 md4 (ref) | 7 | 1 (r122_g1_2 .99/.62) | .154 | union .197 |
- ROUTE offline on recipe16 md4: best leaf F1 .154 -> .151 / .147 / .145; union .197 -> .187 / .185 / .184 (slightly worse). 10.0 core-h.
- Costs OFF 24.0, INH 23.5, SUP 27.5 core-h.

## Continuum SUP (tier 2; md2; ct_route_score.py)
| fixture | base eta2 / W1 | INH eta2 / W1 | SUP eta2 / W1 | SUP leaves (member t) |
|---|---|---|---|---|
| ctpk | .185 / .101 | .255 / .075 | .027 / .125 | 6 leaves .55-.66 |
| ctun | .385 / .133 | .448 / .097 | .176 / .160 | 4 leaves .38 .41 .58 .66 |
- SUP destroys the continuum axis (the swinging domain is outside the order support): splits follow something else
  (z 4-15, all EM-routed). ctun SUP r1: guard SKIPPED. Cost 33.9 / 28.8 core-h.

## recipe18d (tier 3 build) = recipe18c + locality_guard_children_r17d.diff (d6e45d9f): build_r18d, fast gate 15/15, simple_exec f6fb0afe

## r1211 (tier 1; md4, kway 2; leaf_scorecard; baseline ancestral_dsnr r1211_s1_off, no corner leaf in any arm)
| arm | leaves | best closed-leaf F1 (leaf, rows, ref, arm) | notes |
|---|---|---|---|
| base | 5 | .642 (r12_g1_2, 7776, 5297, 0.92/0.61) | |
| INH+PG | 4 | .612 (r12_g1_2, 7906, 5087, 0.90/0.73) | gate 0.09 %; 6 inherited seedings; r11 not split further |
| SUP | 5 | .279 (r21_g1_2, 6289, 2096, 0.51/0.69) | all maps closed arm 0.48-0.56: ARM CUT |
- ROUTE offline on base: best leaf F1 .642 -> .583 (down) / .596 / .595 (flat1/2): worse. 6.0 core-h.
- Costs INH 26.4, SUP 28.1 core-h.

## qtoy2 SUP (tier 2; md2; disc_score.py; base / INH from the continuum session)
| arm | leaves | NMI | best-leaf F1 mean / min | c18 (closed) | c1 c2 c22 c28 c43 c49 c101 |
|---|---|---|---|---|---|
| base r17d | 11 | .684 | .77 / .55 | .95 | .55 .68 .79 .74 .62 .83 .99 |
| INH (build_inh) | 13 | .683 | .77 / .57 | .95 | .79 .57 .79 .75 .62 .83 .87 |
| SUP | 15 | .532 | .58 / .47 | .47 | .71 .47 .61 .49 .81 .52 .59 |
- SUP halves the closed conformer's F1 and loses .15 NMI; 1 guard skip (r11). 46.5 core-h.

## mm3q SUP and INH+PG (tier 2; md2; disc_score.py; base = recipe17d md4 tree, the md2 baseline comes from recipe18d mm3q OFF)
| arm | leaves | NMI | best-leaf F1 mean / min | combos >= .7 | guard skips |
|---|---|---|---|---|---|
| base r17d md4 | 15 | .640 | .67 / .52 | 1 of 8 | 3 of 14 |
| SUP md2 | 10 | .705 | .75 / .48 | 5 | 4 of 11 (r11, r12, r22, r31) |
| INH+PG md2 (gate 0.02 %) | 12 | .673 | .73 / .61 | 4 | 0 of 9 (5 inherited seedings) |
- Per combo SUP: c111 .78 c112 .68 c121 .85 c122 .66 c211 .48 c212 .88 c221 .78 c222 .92; INH: .70 .92 .61 .74 .64 .91 .68 .62.
- mm3q's moving domains are rigid, textured and well ordered, so the order support keeps them: SUP helps here (unlike RNase P).
  Depth differs from the base (md2 vs md4); the fair md2 base is recipe18d mm3q_OFF_s1. Costs 42.1 / 43.7 core-h.

## (a) Guard fix on recipe18d where the guard used to skip (md2, seed 1)
| run | node | z (crit 3.38) | choose | children L/L_ref | decision |
|---|---|---|---|---|---|
| qtoy2 INH | r11 | 5.21 | EM groups | 0.682 | kept |
| qtoy2 INH | r41 | 4.73 | EM groups | 0.511 | kept (0.011 over the bar) |
| qtoy2 INH | r4 | 3.61 | EM groups | 0.639 | kept |
| mm3q OFF | r11 | 8.99 | EM groups | 0.955 | kept |
| mm3q OFF | r31 | 3.51 | EM groups | 0.834 | kept |
- 0 skipped calls left (qtoy2 INH 9 guarded + 3 children-judged; mm3q 8 + 2).
- qtoy2 INH: the same 13-leaf partition as the unfixed build_inh run (leaf numbering only), NMI .683, c18 .95: no change.
- mm3q OFF md2 (the fair md2 base): 11 leaves, NMI .701, F1 mean .74 / min .56, 4 of 8 >= .7. Against it: SUP md2 .705 / 5 of 8
  (level; 4 unguarded splits on recipe18c), INH+PG md2 .673 / 4 of 8 (slightly worse). The earlier gain vs md4 was depth.
- Costs mm3q OFF 44.3, qtoy2 INH 46.4 core-h.

## (b) r131 seed 2 on recipe18d (guard fix in; no skips here, so = recipe18c)
| arm | leaves | corner leaves | best closed leaf (rows, ref, arm) | F1 |
|---|---|---|---|---|
| OFF s2 | 3 (r2 re-test fails) | 1: r1_s2 0.99/0.51 | r1_s2 13332, 1214 | .148 |
| INH+PG s2 (gate 0.34 %) | 4 | 0 | r1_s1 12193, 1083, 0.92/0.51 | .142 |
| (s1) OFF / INH+PG | 4 / 4 | 1 (0.97/0.54) / 0 | .123 / .147 (INH 0.92/0.52) | |
- Root identical at both seeds (z 16.26); r1 and r2 are pole-less, fflex seeds in OFF, inherited in INH. INH loses the corner at both
  seeds (closed arm 0.92 both times; OFF 0.97 / 0.99): REPRODUCIBLE corner loss. Costs 21.5 / 21.8 core-h.

## VERDICT (recipe18c / 18d, 2026-10-29)
- SUPPORT=order: FAIL. Cuts RNase P's flexible arm (all leaf maps closed arm 0.33-0.56 on r131, r111, r1211); qtoy2 closed F1
  .95 -> .47; continuum axis lost. Nulls pass. Level with base on mm3q md2 (rigid domains); L10_12k gain was an unguarded split.
- INHERIT: FAIL on the RNase P corner at r131 (2 of 2 seeds, 0.97/0.99 -> 0.92); r111 corner kept; r1211 .642 -> .612. Helps the
  continuum (eta2 .185 -> .255, .385 -> .448); neutral on qtoy2 (NMI .684 / .683) and L20; slightly worse on mm3q md2 (.701 -> .673).
- POSEGATE=0.165: inert outside the micelle (root flips 0.02-0.34 % on RNase P, nulls, rare, mm3q); on micelle it fires and equals
  fixposes (micnull PASS, mic3 micelle split kept). Safe to have on; only matters for micelle-like data.
- ROUTE: helps discrete toys only (mm3q flat2 1 -> 6 of 8 combos, qtoy2 per routing session); neutral/worse on RNase P (r1211
  .642 -> .595), rare (.78 -> .55), continuum (eta2 -.02 to -.06); OLE reshuffles half its rows on sub-nat margins.
- Guard fix (children path): every former skip judged (0.511-0.955), all kept, leaves unchanged on qtoy2 INH; 0 skips left.
- Guard-skip frequency: recipe17d trees 0 skips on RNase P / IgM / OLE / qtoy2 / 692k; 3 of 14 on mm3q. Switch arms create them:
  INH qtoy2 3 of 12, SUP 1-4 per tree (rare roots, ctun, qtoy2, mm3q). All were choose -> EM groups.
- Cost about 690 core-h (tier 1 ~258, tier 2 ~287, (a)+(b) 134, builds/scoring ~10).
