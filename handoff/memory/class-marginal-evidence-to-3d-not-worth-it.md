---
name: class-marginal-evidence-to-3d-not-worth-it
description: "2026-10-14 verdict: porting flex2D's class-marginalised evidence (FLEX2D_CLS_MARG, mixture of PPCA) to 3D flex states is not expected to help; the 2D version already hurt on hbK15k"
metadata:
  node_type: memory
  type: project
  originSessionId: 3e8ad5bc-323b-4ac6-9e9c-71b71e79fe55
  modified: 2026-10-14T17:21:04.062Z
---

User asked (2026-10-14) whether flex2D's per-particle model (class mean + class basis, evidence marginalised over z,
pose tangents deflated) could extend to 3D flex. The mechanics port directly: flex_pca's probe_solve_mix already forms
the Woodbury G/b/logdet. The answer given was "probably not useful":
- In 2D it already HURT (hbK15k: base 0.828 -> r2m 0.818, r4m 0.812, deflated 0.813). A learned per-class basis absorbs
  the state difference of mixed classes, so marginalising forgives the wrong label.
- In 3D a per-state B_k is fitted on impure states (RNase P soft-EM S2/S4 is a continuum cut per seed), so it fails the same way.
- Soft-EM already uses full per-state maps, which is the useful part. The 3D bottleneck is per-particle SNR and poses,
  not the form of the likelihood.

**How to apply:** don't re-propose state-marginalised 3D evidence as a separation lever. The only use left is model
comparison (one state + flex vs two discrete states) on the S2/S4 continuum.
Related: [[flex2d-findings-2026-09-21]], [[rnasep-softem-real-data-2026-10-13]], [[flex3d-soft-em-proposal]].
