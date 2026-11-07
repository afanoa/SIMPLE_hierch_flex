---
name: flex3d-soft-em-proposal
description: "2026-10-04 proposal handed to its own session: iterated soft EM over flex 3D states (cross-half, pose-marginalised image-likelihood posteriors -> M-estimation -> repeat); only one-shot M-estimation has helped flex 3D so far"
metadata:
  node_type: memory
  type: project
  originSessionId: 95245f2d-90ec-4676-978e-28a5533c3f54
  modified: 2026-10-05T01:38:36.116Z
---

User 2026-10-04: "Nothing has really improved the flex 3D result so far other than the M-estimation"; still considering a
multi-body flex; cores are tight, so toy tests may run on the Mac (Python prototype first, no Fortran builds unasked).

Proposal (handed off as a separate session, task "Prototype soft-EM state refinement for flex 3D"): (1) per-particle state
posterior from a calibrated image likelihood against the state maps, marginalised over pose (local per-state orientation
search + shifts), computed cross-half; (2) M-estimation: every state reconstructed from all particles weighted by the posterior;
(3) iterate. Existing parts: score_states (flex-weights-refine branch, fixed pose + shift, optional local orientation search),
the per-state weight store, refine3D/reconstruct3D flex_weights=yes with FROZEN weights. Missing: step 3 and pose
marginalisation by default.

**Why flex states, not abinitio3D states:** fixed-pose scoring of abinitio3D multi-state maps was degenerate on 10076 (argmax ==
current label 100 % even cross-half) because those maps are not co-registered; flex states are mu + U*m_k in the consensus frame,
so one pose serves every state. Framing: discrete states = flex's latent model with a GMM prior; closed-form per-state marginal
likelihood from b_i, A_i, plus the image term for variability outside span(U).
Related: [[score-states-likelihood-posteriors-toyk]], [[spliceosome-mest-refine3dauto-2026-10-02]],
[[impsamp-toy-campaign]], [[multibody-flex-means-region-pca]].
