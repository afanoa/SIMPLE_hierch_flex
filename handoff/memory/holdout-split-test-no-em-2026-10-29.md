---
name: holdout-split-test-no-em-2026-10-29
description: held-out split test (fit one half, score the other at 40-20 A): NOT adoptable 2026-10-29; no variant (v1-v5) beats production+floor+guard on 19 nodes; docking, EM refits and best-group noise fits each break it
metadata:
  type: project
---
Held-out split test (holdout_on_r17b_v3.diff, SIMPLE_SPLIT_HOLDOUT=10, default off; replays 2026-10-29, route9/hold).
- v1 DOCKED each group map to the other half's consensus (translation up to 3 px). Inside small focus masks that spoils
  one map, so T = log π.
  - recipe17e gate false negatives: qtoy2 r3 closed conformer; r1211 r11 refits.
  - v3 removes the docking. Numbers below are v3, from 17 nodes.
- r1211 r11 (open substates, re-test z 7.05, rho 0.75) still FAILS in every arm (z 1.5).
  - Its maps differ (dsnr 5.7 at 40-20 A), but particles do not prefer their state.
  - The coordinator decides whether that refusal is acceptable.
- With SIMPLE_SPLIT_HOLDOUT_ITERS=0 (group maps straight from the stage-B board labels, undocked):
  - all 8 null seeds are K*=1 without the locality guard (production splits qtoynull s1/s3);
  - every positive except r11 passes; m33 closed pair z 21 vs production 3.66; qtoy2 r3 z 342;
  - r2121212 z 10.6 vs 0.98, but the 0.03 rho floor refuses it.
- Nulls' stage-B labels score BELOW shuffled labels: own-half selected groups do not generalise.
- With within-half EM (3 iterations), shuffled-label refits rediscover the strong real structure, so the null mean rises
  (qtoy2 3.1, qtoy2 r3 0.87 nats per particle) and qtoy2 r3 is lost.
- Cost per split step: no EM about 4.5x cross_em, state_split wall +18-65 %.
- The shuffle-null sd is tiny, so r2121212's per-particle z is only 1.4.
- VERDICT after 19 nodes (v3 / v5 tables in doc 10.20): not adoptable.
  - The shuffled-label null is miscalibrated on heterogeneous nodes: random partitions gain +0.08/pp at qtoy2 r11, more
    than the real impure split.
  - The absolute gain refuses impure real splits and passes a junk-share split within one conformer (qtoy2 r4).
  - v4 (noise fit per best group) is a selection on noise; it lets three nulls through.
**Why:** this fixes the leak found in [[coarse-retest-band-leaks-in-em-band-2026-10-29]]. Stage B refines with independent
even/odd references (lp set, so no low-res blend), so its labels are valid held-out starts.
**How to apply:** for held-out or permutation tests, never refit EM inside the null (it finds real structure from random
starts). Never dock group maps one by one before scoring particles in their own frames. Doc section 10.20 records the floor decision (does the floor stay next to the holdout?) as open for the coordinator.
- recipe17e GATES (production env + HOLDOUT=10 ITERS=0, 2026-10-29): qtoynull s1 holdout z -15.0 FAIL (floor alone would pass, rho
  0.044), s3 z -27.9 FAIL; both one leaf. Real roots PASS: r1211 z 51.9 (rho .071), qtoy2 z 1719, r111 z 287 (closed 25,528 whole).
- IgM landscape vs truth (landscape_score.py): recipe14 covers 51-66 % of the truth hull at W1 0.060-0.065, recipe12b 6-13 % at
  0.081-0.085 (narrow tree inside the truth).
- FALSE NEGATIVES at deeper real nodes (recipe17e): qtoy2 r3 (closed conformer 18 split, rho 1.33) holdout z -30.9 / pp -16.0 ->
  K* 1; r1211 s1 r11 (open substate, rho 0.75) z 0.03 / pp -38.8 -> K* 1. NOT adopted; production split test + locality guard
  stays. Diagnosis asked: own-label vs mixture-gain scoring, band of the difference.
