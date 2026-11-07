---
name: multibody-flex-design-2026-09-19
description: "2026-09-19 investigation verdict (no code): multi-body flex = flex_pca with a GENERATED basis (6 rigid tangent generators + 1 occupancy column per body); E-step and kernel-weighted delivery unchanged; SIMPLE_COV_COMPOSE already embeds a basis from MRC files with no EM; no signal subtraction in the linear stage, on-the-fly polar-domain subtraction only in a nonlinear polish; regions from rigid-consistency segmentation of the refine3D_states finals; solvent prior (34e69f0c4) is production-only and not the RNA limiter"
metadata:
  type: project
---

Doc: `doc/implementation_notes/flex_rna_symmetry_multibody_investigation_2026_09_19.md` (sections 4-6).

**Model nesting:** learned eigenvolumes (flex_pca) / rigid tangent generators `-w_r(e.grad mu)`, `-w_r((e x (x-c_r)).grad mu)`
(multi-body = E1, [[e1-rigid-tangent-subspace]]) / smooth-field generators (deformation: ENM, splines, LDDMM scalar momentum).
Multi-body is the deformation model with a piecewise-rigid field. Same forward model `y = a T P (mu + U z)`, same E-step
(`embed_latents_with_contrast`, `simple_flex_pca_em_embed.f90:15-28`), same delivery. IgG record: rigid domains > ENM at
matched rank (+0.211 vs +0.144); growth/occupancy term beats warp -> add per-body occupancy columns (DynaMight/3DFlex
cannot model composition). Pose jitter becomes the explicit COMMON block of the body latent (fixes the DEFLATE_POSE dilemma).

**Plug point today:** `SIMPLE_COV_COMPOSE=<dir>` loads `flex_pca_pc00K.mrc` + `flex_pca_probe.txt` and goes straight to
the embedding ("no covariance moment, no EM", `simple_flex_pca_model.f90:355-364`) -- a generated basis can be embedded
with no solver change (it orthonormalises and rescales; a real `basis=generated` key would keep A/rad units).

**Nonlinear polish (only bodies beyond ~lp/2 displacement):** per particle, coordinate descent over bodies, each a local
search in its own polar bank (masked body map, direction neighbourhood + in-plane + centre-induced shift) against the
residual minus the other bodies' current projections = B-1 bank lookups with phase factors inside the search. That is
RELION's on-the-fly subtraction ([[relion-multibody-and-symexpand-facts]]) without its sequential single-body likelihood.
Subtracted stacks are never needed; per-body reconstruction in the body frame = deformed backprojection for free.

**Regions (knob-free):** displacement fields between the refine3D_states finals (patch cross-correlation, no aperture
problem; PCs only give the gradient-parallel component) -> rigid k-means (DynDom/RigidFinder) under envelope + contiguity +
minimum mass; body count by the cross-half gate on the body coordinates; fallback `seg8` density k-means. NU tiers alone
gave 29 fragments in August: local resolution is not rigidity.

**Solvent prior:** 34e69f0c4 = soft ridge `lambda_s(1-w)` in the PRODUCTION regularized solve (per-half Otsu weight,
`simple_pcg_solvent_sidecar.f90`); flex M-step has a hard window + Fourier-shell ridge only (`simple_flex_pca_pcg.f90:418-421,
1204`); flex state delivery on pcg uses `reconstructor_pcg` so `set_solvent_prior` is reachable. Measured not binding on the
RNA ([[rnasep-flex-seeds-vs-refine3d-finals]]); order: state delivery first, then a real-space diagonal in the flex operator
with a dilated envelope. Does not help micelle belts.

**Ranked next tests:** (1) RNA on-axis states from the per-state caches; (2) embedding-only iteration ladder scored on
seed contrast; (3) E1 on the RNA with seg8 vs finals-segmented regions, cross-half gate; (4) cnga1 irrep detector + one-
protomer injection ([[flex-symmetry-two-regimes-not-expansion]]); (5) solvent prior wiring; (6) min_state_frac scaled by
nstates, envelope-based flex init.
