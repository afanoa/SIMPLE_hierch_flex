---
name: tree-as-full-mppca-assessment-2026-11-06
description: "user asked 2026-11-06 if hetero_tree can be a full mixture of PPCA (and whether variational PCA / 'Riemannian multipliers' help): yes formally (rank-0 leaves, hard EM, cross-half accept); full version predicted no partition gain; VB/ARD and a Gaussian-latent likelihood cannot be the split rule"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-06T17:37:39.670Z
---

Assessment given in chat, nothing run:
- The tree is a hierarchical MPPCA with rank-0 leaves: node fflex = one PPCA (the K = 1 limit); site EM / seeds = split proposal;
  stA0 / stage A / stage B = hard EM over (state, pose); split test / merge-up = add / remove a component by cross-half test;
  tree_wiener = hierarchical prior. Every node fits fflex before its split test and the prune deletes it, so a per-leaf PPCA
  readout needs no extra fit.
- A full version (leaf W_k with Woodbury-marginal evidence, soft (state, pose) responsibilities, weighted flex_pca M-step) is predicted
  to give no partition gain: [[klem-refine-late-and-mppca-design-2026-10-16]], [[class-marginal-evidence-to-3d-not-worth-it]],
  [[soft-multistate-fixed-pose-verdict-2026-10-30]], [[final-routing-pass-verdict-2026-10-29]].
- Mean-field VB (variational PCA, VB mixtures of factor analysers) does not cure the noise-lock: its E-step scores the posterior-mean
  (Wiener) map plus a term that does not depend on the image (component size, pose, CTF), which cancels between equal siblings.
  ARD and the free energy keep quality and pose-error axes, because those are real variance.
- Gaussian-latent trap: at d' ~1 along a site, excess kurtosis is -0.08 for two states, -0.05 for an even continuum, 0 for a PPCA.
  "Two PPCAs vs one of rank q+1" therefore calls an even swing discrete at ~1e5 particles per node, and separates real two states
  from it only at ~2e5.
- A tangent-space (SO(3) Laplace) pose marginal is the cheap soft pose marginal, but first order holds only while r*dtheta < d/(2 pi)
  (~1.4 deg at 100 A radius, 15 A); node poses err ~4 deg median ([[ppinfo-lever-verdict-2026-11-02]]).

**Why:** this is the paper's integration gap ([[hetero-tree-paper-plan-2026-11-06]]).
**How to apply:** use MPPCA as the paper's formal model, with cross-half acceptance in place of likelihood or free energy. Don't
re-propose soft or leaf-basis E-steps, or VB, as separation levers without a new map-noise model. Leaf PPCA = readout only.
