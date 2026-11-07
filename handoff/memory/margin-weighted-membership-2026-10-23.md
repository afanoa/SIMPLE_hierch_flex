---
name: margin-weighted-membership-2026-10-23
description: "hetero_tree margin weighting (m8c17, uncommitted; FAILED toy calibration, tau* ~3): q = softmax(E + log pop), E = log-sum-exp over a state's table candidates at T = f*pi*wsqsum/(2*pftsz); likelier-only children, weights multiply down; plus discovery re-cut (m8c16) and the stage-A label-compaction fix"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-23T05:14:11.423Z
---

2026-10-23 night, user (via the coordinator): "yes apply the discovery fix and the margin weighting" as hetero_tree defaults; no
commit until the user says. Worktree sharp-proskuriakova-b68280; batch2 ~/SIMPLE_m8c builds 16 (re-cut) and 17 (weighting + fix),
snapshots ~/SIMPLE_stree/snap/m8c16, m8c17 (both gated).
- Scale is the euclid objective's own likelihood, already in the code as posterior_temperature = pi*wsqsum/(2*pftsz)
  (sigma2 per real component, pftsz samples on a half ring), times f = sum(m^2)/box_crop^2 of the 2D soft mask (~0.45 at stage B).
  Replicate flips are NOT the calibration (same image noise -> lower bound on doubt); they are the check: sum 2q(1-q) >= flips.
  fast-tree's r1 numbers: at T ~ 389-600 mean q 0.83-0.88, q < 0.9 for 38-57 %; holds the flip inequality on every schedule pair.
- Mechanics: prob_tab_neigh fill writes prob_wsqsum_partNN.dat; the prob_align master writes prob_state_evidence.txt (iptcl state
  T E(1..K), nats) every iteration; hetero_tree decide() -> child weights.txt / leaf_<name>_w.txt; node reconstructions publish a
  KERNEL state weight set + m_estimator=flex; state_split / site_rank weigh via state_mergeup%rw; refine3D stages and flex_pca stay
  hard. Structural counts stay on rows; nodes log rows / sum w / effective N.
- Twin r2 re-cut result: twin a 0 -> 3 leaves; vs twin b NMI 0.10 (candidates 4/8 matched). Seeded reruns are bit-identical.
- RESULT (cdisc truth, 4 trees 2026-10-23): q is ~3x OVERCONFIDENT, oracle tau* 2.92 / 3.71 (root A split: mean q 0.92 vs 0.80
  correct; q>0.995 bin 6 % wrong) -> stopped per the agreed rule (no tempering); trees change structure without consistent gain.
  Flip inequality holds (7,454 vs 19,971; RNase P r1 569 vs 20,602) but cannot test calibration. Fix candidates: weights off,
  or real pose marginalisation (local Laplace) in prob_tab. Decision pending with the user.
**Why:** honest membership of ambiguous particles without duplicating them (a soft split needs a weight floor = knob).
**How to apply:** validation runs fulltoy_{h16s1,h16s2,w17s1,w17s2} (stree_wtoy.sbatch, eval wtoy_eval.py); if the oracle tau* on the
toy is far from 1, stop and report, never temper. Related: [[factorial-site-layer-port-2026-10-23]],
[[mac-tar-appledouble-breaks-batch2-build]].
