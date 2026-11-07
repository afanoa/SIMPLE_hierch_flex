---
name: flex-pca-pcg-mstep-identity
description: "A production matrix-free PCG volume solver already exists (rec_backend=pcg) with mask-inside-solve, and solve_coupled_basis_exp is exactly its block-Jacobi preconditioner — warm PCG strictly subsumes the flex M-step with a budget-0 bit-equivalence gate"
metadata: 
  node_type: memory
  type: project
  originSessionId: 77fb61e4-4f59-4f34-a8f5-e43daf901228
  modified: 2026-08-31T18:05:33.691Z
---

`src/main/volume/simple_reconstructor_pcg.f90` (4572 lines, `rec_backend=pcg` in
production reconstruct3D, contract `doc/policies/reconstruct3D_pcg_policy.md`)
already ships: real-space mask INSIDE the solve (u = M w, operator M·A·M —
`set_mask`/`apply_normal` ~:2278-2301, `solve_core` tail), per-shell FSC-tau and
graded real-space quadratic priors matrix-free, sampling-density Jacobi
preconditioner with shell-mean floor, knob-free relative lambda, dx/x <= 1.5e-2
stop, and a Toeplitz/NUFFT-Gram kernel operator (~7x faster/iter, ~3.4% interior
error, not reliably PSD — production caps at 2 iterations, rejects maxits>8).

The key identity (the production code documents the K=1 case in writing): the
flex M-step's per-voxel coupled division `solve_coupled_basis_exp` IS the
block-Jacobi preconditioner of the true coupled normal equations, so today's
M-step = "0 CG iterations + M^-1 b". A warm-started 2-3-step coupled PCG
strictly subsumes it; budget-0 reproduces it bit-for-bit (the regression gate).
The normal operator applies through the direction bank at ndir cost, not N:
A u = sum_d P_d^T [S_d o (P_d u)] using the same per-direction S_d pair tables
the banked adjoint's GEMMs already produce (and that the comlin fan M-step
needs — the two designs share their sufficient-statistics layer).

Why it matters: this is the only route that makes the expanding real-space
prior roadmap (solvent mask, density priors, deflation-as-constraint) EXACT in
the fitted model instead of a post-hoc projection. Composition verdict:
comlin-as-preconditioner rejected (not SPD, redundant with free block-Jacobi);
the tracks compose sequentially — fan M-step for cheap ring-resident probe
iterations, PCG for prior-exact volume iterations (stage-2 polish, final pass).
The one real implementation risk: the coupled adjoint-identity test vs
`project_fplanes_mean_basis`. Also: `cov_packed_cgsolve`'s CPU implementation
(`cg_solve_packed`) did NOT survive the EM split — only the env gate and sizing
remain; `opt/` CG modules are nonlinear GSL minimizers, unusable as a linear
driver. See [[flex-pca-volume-free-comlin-mstep]].

**Correction 2026-09-07:** the line numbers above are from the 4572-line solver at 5aee28fe7. On origin/master (124b277f4) the file is
3284 lines: set_mask :654, NEW set_mask_volume :677 (arbitrary envelope support, dev item 5), apply_normal :1067, apply_normal_kernel
:1199, accumulate_rhs_density :1772, build_kernel :2297, prepare_fused_planes :2567, apply_precond :2686, solve_core :3049. The graded
real-space quadratic priors (Q_s solvent, Q_NU band weights / apply_nu_precision) and the LM pose code are GONE from the solver; the
prior surface is hard support + P_tau + lambda_0 only. Do not cite set_nu_prior/apply_nu_precision as available. Details in
[[flex-pca-real-space-priors-proposal]].
