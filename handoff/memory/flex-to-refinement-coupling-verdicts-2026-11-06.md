---
name: flex-to-refinement-coupling-verdicts-2026-11-06
description: "VERDICT (2026-11-06): every flex -> refinement coupling at frozen poses tested in Python on qtoy2 / r131 / qtoynull / olenull (prior on the state argmin, two-distribution fusion, latent-regularised references, node-level held-out gain) loses to the tree's stage B plus decision-time calibration; numbers, traps, what remains"
metadata:
  node_type: memory
  type: project
  originSessionId: 5482e66b-b4b2-4b69-b534-a7b318bf3360
  modified: 2026-11-07T02:11:36.855Z
---

Scripts `scripts/hetero_tree_proto_latent/f2r_*.py` (cluster copy `math_reworks/latent/`, runs under `runs/<key>/f2r_*`), notes in
the side-session folder `f2r_flex_to_refinement_status.md`. All at frozen stage-B poses, sites basis (flex modes + site factors).
- Prior alpha log w on the state posterior, (beta, alpha) by checkerboard held-out: real information (beats a matched stickiness
  control) but fades as references sharpen (qtoy2 it-1 NMI gain +0.040 at 15k -> +0.018 at 60k; r131 full alpha 0 from it 3).
  Frozen at the between-stage calibration it overrides held-out data (-0.004..-0.025 nats/ptcl) and anchors null splits.
  Alpha = 0 on eight null runs. Flex-mode-only weights (what the tree feeds m_estimator=flex) carry nothing; only + site factors.
- Two-distribution posterior (per-state Gaussians over [latent log-ratios, objective gaps]): worse than the calibrated product
  (qtoy2 -0.044 vs +0.077), the objective's gaps are not Gaussian within a state (objective-only learned model -0.11..-0.25),
  overconfident and negative on nulls with no off-switch.
- Held-out latent gain as a discrete-vs-continuum test: f50 0, cont +0.0006, null +0.0005, f15 (15 % minority) +0.028. Measures
  a noisy minority reference, not continuity.
- Latent-span regularised state references (deviation projected onto the span): stops the drift on r131 full (1.37-1.38 over 4
  its vs O 1.34 -> 1.19) but loses signal where the span is incomplete (qtoy2 0.57-0.59 vs O 0.615, C 0.633) and makes null
  deviations look reproducible (FSC 0.22 -> 0.95; the site factors come from the same split). Raw-FSC weighting restores the
  calibration but returns the drift (r131 15k 1.391 -> 1.295, below O). Shrinkage and stability trade off; nothing beats C.
- Every frozen-pose hard-EM chain drifts below the labels it starts from; the proxy cannot show a gain over production stage B
  (the parent's refine=state_eval harness shows the REAL refinement does not drift at frozen poses: the decay was the proxy's).
- In production (merge_sm build): the marginal state rule (prob_state_marg) moves 2-6 % of labels, r131 enrichment 1.450 ->
  1.466, nulls unchanged: nothing. Posterior temperature by held-out prediction: pose-only beta* 0.10-0.20, joint 0.2-0.3;
  but gate trees at 0.15 lost the closed corner and merged a qtoy2 pair, so 0.3 stands. LATENT SEEDING of stage A (flex +
  site-factor mixture labels replacing site_em's seeds, node replays with the tree's own calls): better seeds (qtoy2 NMI 0.57
  vs 0.31) but the finished stage B ties on both toys (qtoy2 0.650 vs 0.656) and loses the closed-row enrichment on r131 r1
  (1.02 vs 1.25 at both seeds); nulls one child. Stage B converges to the same basin from very different seeds.

**Why:** the user asked for a decision on coupling flex to the refinement's posterior under the hard gate (frozen poses, no
latent recalculation, no Woodbury); under a discrete model the objective is sufficient once references are good, so the latent
can only denoise, and that fades with N.
**How to apply:** do not retry these forms. Calibration belongs at decision time (split/merge), not in the iterations. Untested
and only testable in production: RELION-style sum over poses for the state evidence in assign_rows / gather_state_row. Traps:
`prep_crop` accepts any plane cache of matching shape (a second subsample silently reuses the first); `calib_flexprior2.py`'s
objective-only row is scored at the calibrated chain's references (not an independent chain). Replaying a gate node's steps
needs `SIMPLE_SPLIT_CONTRAST=maps` exported for the split step (`gates/<arm>/split_env.txt`) and `SIMPLE_SEED` as the tree's
(launcher default 1; r131 ran at 2); fixture and child-node projects carry the parent's stacks, so `prg=selection` refuses
them (`write_states` in merge_sm writes states tree-style). Related:
[[flex-frozen-poses-only-no-woodbury-2026-11-06]], [[soft-multistate-fixed-pose-verdict-2026-10-30]], [[inode-prototype-2026-11-06]].

2026-11-08 revisit with sparse-PPCA site posteriors (precision A) on IgM, in `state_eval`:
- A fixed prior helps where the seed carries the arm (10 % +.01-.05 NMI) and anchors bad seeds.
- Its held-out alpha stays > 0 on continuous-arm nulls.
- See [[fixed-latent-prior-on-state-eval-verdict-2026-11-08]].
