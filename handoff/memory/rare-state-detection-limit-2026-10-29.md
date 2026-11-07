---
name: rare-state-detection-limit-2026-10-29
description: "Isolated rare discrete state (qtoy SNR, closed-arm effect) at the tree root: found only at n/sqrt(N) >= 16-17 (1.9-2.1 % at 692k), RE-MEASURED 2026-11-02 under the variant D subspace site test = unchanged; lost at the arm site's endpoint acceptance (r 0.14 vs 0.275 at n_eff 116) then the guard merges the fallback split; each unnecessary binary level costs sqrt(2) in share (verified); lever = climb-by-scoring at the site"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-30T02:29:35.157Z
---

2026-10-29, detection-limit session (rare_state_detection_limit_2026_10_29.md + rare_state_detection_curve_2026_10_29.png +
hier_ntests_on_r17d.diff md5 482267a5 in the shared notes; scripts hetero_tree_20261028/detlim/scripts; 25 runs, ~330 core-h).
- conf18 (rel 0.60): FOUND 20 % of 6k (n/sqrt N 15.5, F1 .78) and 10 % of 60k (24.5, found at depth 1, F1 .73); root bar 10 % of
  12k (11) and 5 % of 60k (12.2); discovery loss at <= 7.3 (0.25-3 % at 60k, 2/4 % at 30k): flex latent blind, no site pole, no
  route A candidate. conf22 (0.43) lost at 15.5; conf02 (0.21) lost at 24.5. Count and fraction both matter: n/sqrt(N).
- minpop, choose and merge-up never decided a loss. Found cases came via the fflex fallback + enrichment, never route A.
- SIMPLE_SPLIT_HIER (node-size Bonferroni, default off): passes the bar but the guard merges the diluted split (0.397/0.487);
  quality splits <= 0.34, so a guard tweak would be a 0.03-0.06 margin. All 4 nulls stay 1 leaf under HIER.
- Correlated axis (conf18 beyond the end of a qtoy2-like swing majority, 2 / 1 / 0.5 % of 60k, maxdepth 2): lost in all three.
  The swing branch enriches the share (2 -> 2.4 -> 6.5 -> 12.2 %) but n/sqrt(N) only goes 4.9 -> 8.5. Its companions are the
  swing neighbours (conf18 vs conf22 rel 0.41 = the medium effect). A200 is lost at the bar at depth 2 (z 2.38); the deepened
  child is lost at discovery. The isolated rule holds inside the branch, at the effect against the branch's own majority.
- Trap: finished trees prune split/s.simple, site masks and stage-B maps; replayable runs need SIMPLE_TREE_KEEPMAPS=yes.
- detlim/ holds ~147 GB; 14 GB of noise-free stacks moved to detlim/_to_delete_20261029 for the user.

**Why:** the user asked (2026-10-29) to investigate the detection limit.
**How to apply:** RNase P closed (1.26 %) is found only because it lies on a shared motion and is enriched down the tree; an isolated
state below n/sqrt(N) ~15 will not be found. Related: [[rnasep-closed-state-must-be-retained]], [[null-split-diffuse-real-split-localized-2026-10-28]].

**2026-11-02 re-measurement (session tree-rare-state-limit; note doc/implementation_notes/tree_rare_state_limit_2026_11.md on branch
claude/tree-rare-state-limit; batch2 detlim2/, 26 roots, 218 core-h, 235 GB kept maps):**
- Arms old = 9005a8c77 (halves test, dist_sdisc/build12) vs new = claude/tree-dist 26c350abb (variant D subspace test, detlim2/build_vd);
  pool2 fixtures 3 / 4 / 5 / 6.5 / 8 / 10 % of 60k + null, seed 2 at 6.5 / 8 %, random halves of the 8 / 10 % roots (30k); recipe18d keys,
  minpop 260 (130 at 30k), maxdepth 0, SIMPLE_SPLIT_CONTRAST=maps, KEEPMAPS. Both arms IDENTICAL as delivered leaves on all 26 runs.
- Limit: lost at n/sqrt(N) <= 15.9, found at >= 17.1 -> c = 16-17 = 1.9-2.1 % at 692k (13.2-14.3k particles). Unchanged by the new test.
- Losing stage: discovery finds the arm site from 4 % on (0.33-0.60 of the rare-difference energy); the site EM tracks the state
  (posterior AUC 0.75-0.89) but the endpoint's cross-half deviation fails acceptance (r 0.13-0.15 vs 0.275 = 3 Fisher SE at the site's
  n_eff 116; halves test 2.1 SE) -> no pole; fflex fallback enriches to purity 0.14 / recall 0.90; re-test PASSES (z 5.9-6.2) and the
  locality guard merges (L/L_ref 0.33-0.48 < 0.5; 0.60-0.86 above the limit; nulls 0.21-0.24). Guard, not the bar, sets the limit now.
- Arms differ only at the pole step at 8-10 %: old seeded twice from a non-rare pole (leaf purity 0.19-0.25); new formed a 717-seed pole
  from accepted arm endpoints at 8 % (seed assignment: a 10,000-N* class drew 59,145 seeds -> flagged to the site_rank owner) and used
  the fallback; new leaves purer (0.30-0.44); new 15 % cheaper per root.
- Oversplit cost sqrt(2) per unnecessary level VERIFIED: 8 % found at 19.6, both halves lost at 13.7 / 14.0 (bar; guard 0.482); 10 % found
  at 24.5, both halves found at 17.1 / 17.5. At 692k: 1.9-2.1 % root -> 2.7-2.9 % one level down -> 3.8-4.1 % two down.
- Lever (proposed, not implemented): climb-by-scoring at the arm site (one more site_em scoring pass from the endpoint's own posterior
  before the subspace test; audit 6.5) to raise endpoint purity ~2x -> r above 0.275. Do NOT lower the guard (0.1 margin to the nulls)
  or touch the re-test. Oversplit criterion proposed: guard + P/N >= 2 floor (merges qtoy2 r4 only; closed corner margin 1.35x).
