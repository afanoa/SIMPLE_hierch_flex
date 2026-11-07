---
name: igm-arm-flex-prior-verdict-2026-11-07
description: "VERDICT 2026-11-07 (guide 7.7, IgM fixture arms a1/a2 + null twin): the focused flex coordinate's per-particle posterior as a refine=state_eval prior adds no arm information beyond references built from the same coordinate (+0.003-0.006 NMI at best, costs accuracy; with oracle references it hurts a2 0.144 -> 0.116); oracle references at the root's poses are the headroom"
metadata:
  node_type: memory
  type: project
  originSessionId: 77c8a7fa-5cd0-4244-b81c-f308ec6cd639
  modified: 2026-11-07T07:57:34.247Z
---

Handoff: `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/igm-prior/VERDICT.md`. Runs: `math_reworks/latent/igm_prior/{fix_a1s,fix_a2,null_a1s,null_a2s}`.
Scripts: `latent/igp_*`.
Setup: the root's stage-B call at state_eval, pcg_mskfile = the true arm mask, root poses (n.simple). References are rebuilt
through that call from fixed labels (latent split or truth). The prior is the cross-half GMM log posterior of T1's best focused
axis.

- **a1 (K 2).** Latent 0.144 NMI / 0.639 accuracy. Latent-built references, no prior: 0.168 / 0.714. With the prior: alpha 1
  0.171 / 0.687, alpha 2 0.165 / 0.670. Oracle references 0.197 / 0.755; with the prior 0.203 / 0.734.
- **a2 (K 3).** Latent 0.105 / 0.471. References, no prior: 0.103 / 0.473; with the prior 0.109 / 0.475. Oracle 0.144 / 0.612;
  with the prior 0.116 / 0.534. The prior drags labels to the latent's even-thirds split.
- **Null.** The prior overrides the arm objective on 98 % of particles and removes, rather than adds, structure: a1 vs angle
  bins 0.0046 -> 0.0002; a2 0.054 -> 0.051 / 0.048.

**Why:** the coordinator's guide 7.7 asked whether the latent helps where per-particle evidence is weak; the answer is no
here as well as on the gate fixtures.
**How to apply:** do not propose per-particle latent priors for arm states again. Purer references (+0.03-0.04 NMI over
latent-built ones at fixed poses) and arm-shaped masks are the route. Trap: my first 1-D EM, initialised at quantiles,
split a1 on a tail (NMI 0.0004). Use T1's sklearn GaussianMixture (n_init 3, seed 1) to reproduce its axes. Related:
[[stageb-oracle-references-limit-2026-11-07]], [[flex-to-refinement-coupling-verdicts-2026-11-06]], [[site-mask-rework-direction-2026-11-07]].
