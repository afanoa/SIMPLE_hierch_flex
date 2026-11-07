---
name: split-test-holdout-and-drain-traps-2026-10-21
description: "split test (2026-10-21): hold-out re-test beyond the EM band, docking, pairwise equal N, joint decision at consensus + node-refined poses, SSNR junk; traps (drain bug, imprint, shift-error sorting, merge power); final rflex6split run rflex_s30split5 12 leaves NMI 0.45; hold-out blind to low-res states (f692u7 closed arm)"
metadata:
  node_type: memory
  type: project
  originSessionId: eaa5de8d-74b0-4feb-8676-4b043e97e641
  modified: 2026-10-23T23:07:56.795Z
---

Found on rflex_s30a nodes + panel_null_50k / panel_two_50k (rnasep_klem/smdiag, routes/toy55 smdiag*.sbatch, splitbench10/11.py):
- Candidate states derived from the data (flex -> refine3D) carry a selection imprint shared by both halves; the cross-half EM
  SUSTAINS it in the EM band (r11: deviation FSC grows 0.8 -> 0.9 at 58 A over 3 iterations, NMI(conformer) stays 0.008). So the
  post-EM re-test over 39-10 A is inflated (r11 0.42-0.88 vs null 0.14); only `--retest_beyond_em 1` (labels from <= 20 A, test
  16.7-10.2 A) is valid. State-specific refined poses are NOT the main cause: r11 at consensus poses still 0.46-0.76 before the EM.
- DRAIN BUG (splitbench9/10 cross_em): a group drained in the last iteration keeps its MAP particles -> re-test equal N = 77 ->
  no power, everything merges (two_50k hold-out K* 1; the earlier null_50k "3 = truth" with the hold-out was this artifact).
  splitbench11 returns argmax over alive groups. Also em_min 1000.
- My splitbench11 first version took the merge-loop statistic from the LAST pair_and_test call (the full-band log call):
  every hold-out decision used the full band. Fixed (pst captured after the band call).
- Amplitude spectrum matching to the consensus amplifies noise where a group has no signal: blob junk merges into the conformer.
  Blob junk IS a blurred consensus, so any blur-invariant test merges it: screen junk BEFORE the test.
- Junk screen metrics (cross-half, vs the node consensus re-referenced without flagged units), bands 39-19.5 / 19.5-10 A:
  projection blob 0.20-0.31, broken (half molecule) 0.47-0.69, real units 0.94-1.07; per-particle SSNR blob 0.08-0.16, broken
  0.5-0.6, real 0.74-1.3. Rule projection < 0.75 or SSNR < 0.35: null_50k flags exactly the 4 junk units (10.4k of 10k junk).
- Pose error by group (pose_err_groups.py, truth all.txt): r11 groups do NOT differ in consensus pose error (means 0.01-0.05 deg);
  the r11 refined poses are worse than the consensus (2.7 vs 1.4 deg median) because the sim's consensus = truth + small error.
- r11 at CONSENSUS poses kept K* 4-6 even with the hold-out: the cross-half EM sorted particles by the 3-D direction of their
  consensus SHIFT errors (pose_err_groups.py, 3-D frame: group means differ by up to ~1 px, se 0.01; rotation errors equal) ->
  group maps = translated copies, reproducible at all shells (no cliff at the EM-band edge). Truth composition (group_truth_diff.py)
  predicts only 4e-5..2e-3 relative difference energy (r12: 0.06-0.16). FIX: splitbench12 --dock 1 (3-D translation of every group
  map to its half consensus over the EM band, in the EM references and the test): r11 -> K* 1 (equal N; full N kept 2 at 0.05-0.07 =
  leakage decaying from the EM-band edge -> --rb_margin 4), null_50k 1 + junk, two_50k 2 + junk. Not lineage (NMI parent 0.007).
- FINAL test (splitbench13, 2026-10-21): junk screen (proj 0.75, SSNR 0.35) -> cross-half EM <= 20 A (refs docked on the EM band)
  -> re-test shells kb+2..Nyquist(crop 48) with every group map docked by a 3-D translation fitted ON THE TESTED BAND (--dock 2;
  a shift fitted on the EM band is biased by the imprint the same way in both halves) -> max of 20 stratified shuffles; equal N.
  Results: null_50k 1 + 4 junk units, two_50k 2 + junk (NMI 0.73), r11 1 (was 8), r12 4 (holo/apo families). Margin 4 is NOT more
  conservative (narrower band -> lower null max -> r11 2). Recursion: routes/toy55/rflex6split.sbatch + rflex6split_decide.py
  (rflex5.sbatch belongs to another session); rerun rnasep_klem/rflex_s30split2 (job 54039392, REUSE root of rflex_s30a).
- FINAL (2026-10-21 15:09): splitbench16 = splitbench15 + pairwise equal N (every pair at the N of its smaller member: a small group
  vs a big merged one picked up the residual shift nuisance, 0.05 FSC); docking to the OTHER half's consensus (own half shrinks the
  shift), merged groups re-docked; junk = SSNR < 0.35 only, low projection + signal = 'partial' (real RNase P: 4 compositional units
  were flagged by projection). Consensus poses keep a shift-sorting nuisance (r11 2), node-refined poses a reference bias (panel mixed
  unit, null conformer split) -> rflex6split runs BOTH and keeps a split only if both support it (union-find in rflex6split_decide.py):
  null_50k 1 + partial + junk, two_50k 2 + partial + junk, r11 1, r12 7, r21 3, r22 3 (swing kept). Rerun rnasep_klem/rflex_s30split5
  (reuse links rflex_reuse6; nodes of rflex_s30split2(_rN)): 12 leaves (1 junk, 1 partial), NMI 0.45 (old 34 leaves 0.42), purity 0.36
  (0.41), swing eta2 0.22 (0.38), extent 10.1-17.9 A in 2 states (8.8-22.4 in 7), blob junk precision 0.994 recall 0.51 (broken =
  partial). OPEN: the band hold-out is blind to low-res/localized states (real f692u7 closed arm 0.87 merged; smdiag12 tests --focus
  and --em_lim 30). Variant runs: rflex_s30split3 (equal-N merges, 10 leaves NMI 0.45, swing 2 states).
- SUPERSEDED 2026-10-23 (m8c39, e22d3d584): the em_min 1000 floor is GONE from the Fortran state_split. Once the re-test used
  per-pair equal N, the floor was only a size floor. It dropped real minorities in 42 tree splits; in 17 of them every group was
  dropped, giving K* = 0. CAHRA r2: 998 / 986 dropped -> K* 1, now 2 on both seeds. Now a group leaves only when one half holds
  none of it. 31-node replay vs truth: 7 added toy splits supported (junk, an IgM arm), 1 not (p27ct1 r212, statistic at the
  null max), null50 unchanged. Never reintroduce a size floor (the user rejects them).
**Why:** a split test on data-derived labels needs labels and test band disjoint in frequency, and every reported K* needs the
EM's drained groups gone. **How to apply:** hold-out re-test beyond the EM band + junk screen, NO em_min floor (see above);
check logs for "drained" and the re-test equal N. Related: [[klem-split-test-statistic-traps-2026-10-19]],
[[klem-cross-half-final-pass-2026-10-19]], [[klem-null-leaves-sort-by-alignment-quality-2026-10-19]].
