---
name: flex-symmetry-two-regimes-not-expansion
description: "flex_pca symmetry: pgrp=cN REPLICATES every particle at all N mates in every accumulate (symmetric basis enforced, symmetry-breaking impossible); pgrp=c1 (standalone default) fits a C1 basis on FOLDED poses (asymmetric modes = N orbit copies at 1/N variance, register scrambled); refine3D_states flex=yes copies pgrp so the whole Cn multi-state pipeline is symmetric-only; neither is symmetry expansion"
metadata:
  type: project
---

Code: replication loops `simple_flex_pca_pcg.f90:618-660` (pcg kernels), `simple_flex_reconstructor_latent_ops.f90:104-160`
(gridding coupled insert), `simple_flex_gpu.f90:428-471`; polar bank gated to nsym==1 (`simple_flex_pca_em_iter.f90:288-303`);
`pgrp` cmdline-only, default c1 (`simple_parameters.f90:334`); `refine3D_states` does `cline_flex = cline`
(`simple_commanders_refine3D.f90:1223`) so pgrp reaches the flex init. Refine3D_states then searches symmetric
references at poses folded into the asymmetric unit -> protomer-specific states cannot be found by construction.

**Why it matters:** "cnga1 (c4) produces little heterogeneity" is consistent with this, but cnga1 is the wrong
test (single apo/closed population, belt variability -- [[cnga1-is-the-apo-closed-dataset]]).

**Design for symmetry-breaking modes (not implemented):** (1) detector on the C1 basis: irrep energy fraction per
component + cross-half reproducibility + orbit-closure statistic (seconds on box-64 eigenvolumes); (2) register-
marginalised E-step: evaluate the N poses R_g R_i per particle against a C1 basis, responsibilities from the
heterogeneity residual only (the symmetric mean term is identical for all g), M-step accumulates at each register
with its responsibility; flat responsibilities reproduce today's replication, sharp ones give a register-resolved
C1 fit; the symmetric solution is a fixed point so seed asymmetrically (one orbit copy, others deflated); cost N x
E-step, no particle duplication; (3) for refine3D_states: RELION's route = particle_symmetry_expand (rows copied
verbatim incl. random subset, only Euler angles recomposed) + C1 refinement + one-protomer focus mask; SIMPLE has
neither program nor a public focus mask; (4) validate ONLY with an injected one-protomer displacement on real
images ([[pfcrt-injection-floor]] machinery). Power is bounded by per-particle register discriminability at the
covariance band: bound factors / one-protomer gates are realistic, subtle protomer tilts are not.
Related: [[cnga1-c4-pose-folding-scrambles-protomers]], [[relion-multibody-and-symexpand-facts]].
