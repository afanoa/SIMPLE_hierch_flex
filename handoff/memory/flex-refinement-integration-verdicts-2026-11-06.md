---
name: flex-refinement-integration-verdicts-2026-11-06
description: "VERDICTS 2026-11-06: frozen-pose flex <-> refinement coupling tested both ways on qtoy2/r131/qtoynull/olenull; nothing worth building; site factors are the only refine -> flex signal"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-06T21:21:53.002Z
---

Tested at frozen stage-B poses for the 3D tree paper's "integrate the latent model with refinement" question
[[hetero-tree-paper-plan-2026-11-06]], inside the gate [[flex-frozen-poses-only-no-woodbury-2026-11-06]]. Python prototypes in
`scripts/hetero_tree_proto_latent/` (untracked), runs in `math_reworks/latent/runs/fast_<key>/` on batch2; 15k subsamples plus
full-node replications (the f2r session's tables: `runs/<key>/f2r_<tag>/f2r.log`).

- flex -> refine (log-prior alpha*log w_flex on the refinement's state choice, (beta, alpha) by checkerboard held-out likelihood):
  gains only with site factors in the latent (qtoy2 15k NMI +0.04, replicated on seed 2); the full 60k node halves it; r131 full
  node calibration switches it off (alpha 0 from it 3); frozen at its it-1 calibration it holds closed enrichment 1.42 vs 1.375
  but overrides held-out data (-0.025 nats/ptcl) and keeps the olenull split. Not built.
- refine -> flex: refitting flex's class Gaussians from the refinement's posterior each iteration gives no gain and turns olenull
  alpha 0 -> 1 by it 4 (self-reinforcing), rejected; objective reliability weights are ~1 everywhere and anti-correlate with SNR
  (qtoy2 r -0.24), dropped; the real flex_pca rerun at stage-B poses lifts its own-latent ceiling (qtoy2 NMI 0.477 -> 0.504) but
  changes nothing once site factors are in or downstream, nulls unchanged, not wired into the tree.
- Site factors (recB0 state deviations inside site masks) are what lets the latent see the refinement's states: qtoy2 ceiling
  NMI 0.48 -> 0.62, r131 closed AUROC 0.65 -> 0.69.

**Why:** these couplings will be proposed again for the paper; they cost a day to test.
**How to apply:** a new coupling must beat these numbers on the same four fixtures first. Never use a Python PPCA basis refit as a
flex proxy (it fits noise: latent prior scale ~1e31 on r131); rerun the tree's own fflex call. Traps: flex_pca refuses a
multi-state project (one population); make one with `simple_exec prg=selection oritype=ptcl3D infile=<one 0/1 per row> mkdir=no`;
a node_proj=yes node project fails selection's STK/PTCL count check; `prep_crop` reuses any plane cache whose shape matches.

**Decisive, in the real refinement (2026-11-06 night):** `refine=state_eval` (branch `claude/state-eval`, commits 8fb58408e, 3018c2808,
eb58a5455; build `hetero_tree_20261028/state_eval/install`) reassigns states at the stored poses with stage B's own call and maps; it
reproduces stage B (qtoy2 99 % of labels, NMI 0.656 -> 0.653, no drift over 4 its; the Python stand-ins' decay was their own). With the
flex + site-factor prior (`SIMPLE_STATE_PRIOR`, alpha 1 / 2; scripts `latent/seval_*`): qtoy2 NMI -0.001..-0.009, r131 closed enrichment
-0.01..-0.02 (r1) / flat (r2), toy holo F1 +0.02..0.03 and closed unchanged, nulls lose 3-12 % of labels. Not worth building.
Gates at prob_beta 0.15 (pose-only held-out optimum; master draws poses only, the state is the argmin): r131 closed corner lost
(0.91/0.53), qtoy2 NMI 0.748 with the 1|2 pair merged; 0.3 kept it (seed 2), 1.0 lost it.
Latent seeding of stage A (flex + site-factor latent labels instead of site_em, the node replayed with its own calls; f2r session):
seeds start better (qtoy2 NMI 0.571 vs 0.309) but the finished stage B converges to the same basin (qtoy2 0.649 vs 0.656, toy 0.172
vs 0.168) or worse (r131 r1 closed enrichment 1.02 vs 1.25); nulls one child. Stage B is nearly seed-invariant. Not worth building.
