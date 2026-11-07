---
name: flex-pca-volume-free-comlin-mstep
description: "The volume-free comlin M-step proposal (doc/for_developers/ideas/flex_pca_volume_free_mstep_proposal.md) — final architecture after the 2026-08-31 adversarial design pass, with the refuted claims listed"
metadata: 
  node_type: memory
  type: project
  originSessionId: 77fb61e4-4f59-4f34-a8f5-e43daf901228
  modified: 2026-08-31T18:24:24.626Z
---

Proposal doc: `doc/for_developers/ideas/flex_pca_volume_free_mstep_proposal.md`
(2026-08-31), based on Van et al. IUCrJ 2026, 13, 354-363 (eq 13/15 comlin
reprojection, eq 19-20 ring FRC). Deleted single-channel implementation of the
paper's math is at git `b818f7000^` in `src/main/pftc/simple_polarft_ops_restore.f90`
(+`_geom`); the "common line PPCA design doc" the user remembered is
`fourier_domain_latent_eigenvolumes.md`, recoverable at `5219186fd^`.

Final architecture (NOT the naive volume-free loop): per-direction C_d/S_d tables
(the banked adjoint's GEMM outputs, kept instead of splatted) -> comlin fan
accumulation WITH the self term -> per-ring-coefficient KxK coupled solve ->
per-iteration BAND-LIMITED SMALL-VOLUME realization (40-64 cubed, ~1-3 s) where
mask3D_soft/deflation/orthonorm run exactly as today and which feeds the exact-
direction E-step; full volumes only at stage boundaries. Honest wins: serial
per-voxel solve + FFT tail + exp grids + bank rebuild leave the loop (~2-2.6x
probe, 1.4-1.6x single run); the multiplicative case is the halfset ladder
(3-6x) and GPU residency.

REFUTED in the design pass (do not re-propose):
- "solved rings are the next E-step bank" — direction quantization is systematic,
  not noise; ndir^(-1/2) decay, ndir also inflates cost quadratically; realize
  small volumes instead (BANKFWD ARI 0.836 precedent).
- "ring FRC replaces the volume FSC as SSNR" — under a shared basis the E[z*n]
  self-interaction is common-mode across e/o halves (Einstein-from-noise); valid
  only between independently fitted halves -> the halfset-EM sibling is a
  precondition. LODO/half-set CV give ZERO solvent protection (coherent solvent
  reproduces).
- boundary-only masking / the crude two-stage — solvent takeover 1/(1+c*rho^2J)
  plus poisoned K* on the ladder; in-loop mask is load-bearing (projected EM).
- eq 15 leave-self-out as an M-step — it is a cross-validated reprojector; the
  self-plane term must be ADDED back or the fixed point is a smoother.
- the deleted code's DOUBLE radial Jacobian — correct only where thin-line
  num/den weights cancel; band gathers need a single source-side k.
- comlin-as-PCG-preconditioner — a fan sweep is not SPD and is redundant with
  the free block-Jacobi.
- comlin operator applied per-particle in the E-STEP — it is the SAME linear map
  as realize+sample (P_i o fan-backprojection identity, value-transport form),
  7-15x dearer, worse kernel, and mask-incompatible; the realization IS the
  operator, amortized (proposal SS4.3).
- LODO E-step as latent de-biasing — removes ~O(1/ndir) of the shared-U FRC
  inflation (the common-mode field is direction-smooth and re-synthesized from
  neighbors); LODO's honest role is the eq-19 FRC construction only (SS3.3).
- per-direction (anisotropic) tau^2 shrinkage — sigma_FRC ~0.14-0.2 at fan
  per-(direction,shell) counts makes it order-unity noise; aggregate per shell.
Adopted from the paper in SS3.3: eq 19-20 shrinkage with DIAGONAL tau^2_q scaled
against shell-sum A_qq (not bare ctf^2), FRC on LODO rings; v1 keeps post-solve
Wiener (fixed-point-preserving), in-solve tau^2 is its own A/B convergent with
PCG Phase C and sibling SS7 cross-fit shrinkage.
Cost model correction: O(ndir^2 * nk * K^2) (pair tables don't factor), fan pair
state ~3.7 GB sp at K=8 (K^2 wall at K~16-24); fan has a 2-5x per-coefficient
variance penalty vs the KB splat.

Phase 0 is pre-registered in the doc: P0a mask-load-bearing arms run on the
CURRENT binary (~60 lines of env gates, RNA mskdiam=210 + IgG + 10076) BEFORE
any comlin code; P0b operator-replacement harness (K=1 FRC>=0.98 arbitrates the
Jacobian); P0e z-permuted FRC null. Obsfield nearest-cell field failed its own
historical quality validation (recorded in polar_observation_backed_restore.md
at b818f7000^) — optional P0b arm only. Related: [[flex-pca-pcg-mstep-identity]],
[[flex-pca-halfset-em-varimax-proposal-docs]], [[flex-pca-em-speedup-plan]],
[[flex-pca-solvent-mask-dominates-rna-basis]].
