---
name: solvar-marginalizes-z-during-fitting
description: "SOLVAR's ML objective IS the z-marginalized marginal likelihood, so the \"SOLVAR does not integrate out the conformation\" novelty claim in flex_polar_joint_pose_heterogeneity.md §3.4 is wrong"
metadata: 
  node_type: memory
  type: project
  originSessionId: 7c9ef785-d2dc-4810-888d-3b610822152a
  modified: 2026-08-06T17:27:41.334Z
---

SOLVAR (Yadgar, Lederman & Shkolnisky, arXiv:2602.17603, 20 Feb 2026,
code at https://github.com/RoeyYadgar/SOLVAR) never estimates `z_i` during fitting.
It starts from the marginal image law `Y_i ~ N(P_i mu, P_i V V* P_i* + sigma^2 I)`
and reduces it with Woodbury + the matrix determinant lemma to
`M_i = I + (1/sigma^2)(P_i V)* (P_i V)` (r x r), giving quadratic term
`(1/sigma^2)[||y||^2 - (1/sigma^2) y*(P_iV) M_i^-1 (P_iV)* y]` and log-det `log|M_i|`.
`z_i` is recovered post-hoc by a closed-form ridge solve after V and poses are frozen.

**Why:** [flex_polar_joint_pose_heterogeneity.md](doc/for_developers/ideas/flex_polar_joint_pose_heterogeneity.md)
§3.4 (line 181) claims "the conformation can be integrated out — this is the part
SOLVAR does not do" and calls it "strictly stronger than what SOLVAR does" (line 21).
That is the same objective SOLVAR's Eq. 9 already uses. The repo's own
[recovar_improvement_directions.md](doc/for_developers/ideas/recovar_improvement_directions.md)
line 631 records `M_i` correctly, so the two docs contradict each other.

**How to apply:** Do not claim z-marginalization as novelty over SOLVAR. The genuine
differentiators left are: discrete polar direction search + one-FFT in-plane instead of
autodiff gradient descent on rotation vectors (SOLVAR needs differentiable interpolation
— it had to drop nearest-neighbour slicing for trilinear 2x-oversampled, and says NN
makes pose refinement impossible); and independent even/odd halves through the pose loop,
which SOLVAR has no equivalent of. See [[flex-pca-latent-is-997-percent-noise]] —
SOLVAR's design implicitly concedes that point estimates of z are untrustworthy.
