---
name: flex-pca-real-space-priors-proposal
description: "2026-09-07 proposal for real-space priors in flex_pca via PCG (doc/for_developers/ideas/flex_pca_real_space_priors_pcg_proposal.md): Phase 0 = molecular envelope replacing the mskdiam sphere at seven tail/delivery sites, no solver code; origin/master already has set_mask_volume/pcg_mskfile; merging it into the flex branch conflicts in six files (sigma2 refactor)"
metadata:
  type: project
---

Proposal doc (untracked on the Mac master worktree): `doc/for_developers/ideas/flex_pca_real_space_priors_pcg_proposal.md`.
Supersedes §5.4/P0a/P0d of the volume-free M-step proposal; absorbs brainstorm §4a/4c; retracts 4e (Laplacian, removed in production).

**What the latest code already has (origin/master 124b277f4, 2026-09-07):** `reconstructor_pcg%set_mask_volume` (arbitrary
[0,1] envelope as the solve support, dev item 5, ca450fe54), the pcg-gated `pcg_mskfile=` override, `set_pcg_solve_support` and
`build_pcg_state_support` (conservative density envelope = `automask3D` at ENVMSKLP_DEFAULT 20 A, ENVMSKWIDTH 7 voxels) in
`simple_rec3D_pcg_strategy.f90`, base-solve warm starts + cold restart on indefiniteness. REMOVED there: Q_NU replay precision and its
auto-lambda/auto-target controllers (2026-09-06, "did not generalize"), the graded NU band-weight operator, and the LM pose code
(solver file 4572 -> 3284 lines). So the production real-space prior surface today is HARD support only; soft (+lambda_m(1-M)^2)
is documented-deferred because a strength knob is the auto-lambda failure mode.

**Branch geometry:** `box/flex-pca-conventions` (cbcd762bd, E2; box checkout a21fac77e G2) and origin/master diverged at dc2eb9a74
(62 vs 46 commits). `git merge-tree` dry run: 6 conflicts — simple_sigma2_files, commanders_flex_pca, em_fit, em_iter, em_pose, model —
all from upstream's canonical-sigma2 refactor (0d7b5f09a, 95d4adbfd) hitting files the alignment steps moved; PCG solver/strategy
merge clean. `box/h0901-int2` holds 5 commits not on conventions incl. d9eb70d34 `SIMPLE_COV_FOCUS_MSK` (arbitrary post-solve mask on
merged basis + both stashed halves) = the Phase 0 hook to extend. Plan: branch `flex-pca-realspace-priors` from conventions, merge
origin/master, cherry-pick d9eb70d34+f3d0d9c92, regen args (lib64 trap), rebuild, re-baseline gate_base.

**The seven sphere sites (box branch):** em_mstep.f90:287/319 (tail), :507 BG template, :544/:584 deflation templates,
em_pairmerge.f90:401, em_basis.f90:88, em_fit.f90:575, rec3D.f90:510-522 + 568-569 (delivery + eo-FSC copies). The mean is never
masked; the 2D image disc (COV_MASK_IMAGES) stays off. Phase 0 = `SIMPLE_COV_ENVMSK=0|1|<file>`, envelope from vol1 via automask3D
(Otsu non-tight, largest CC, dilation scaled to box_crop, cos edge) then `E * sphere(msk)`; built at vol1's box and Fourier-resampled;
vol1 smaller than box_rec -> sphere at delivery. Arms A0-A5, positive control = RNA mskdiam=210 (PC1/PC2 concentration 0.02 -> lobed),
10028 20k harness d'/purity + R2(k=5) >= 0.208 without a smaller sphere. Decision rule: default flip only if no dataset regresses beyond
its rerun floor AND the positive control flips; H0 lowers Phase A below Phase C (state maps on the production PCG solver, needs the
pending "weight both B and D" hook).

**Phase A design:** pair Gram kernel Khat_qr accumulated at DOUBLED coordinates in the coupled insert (same size as rho_e; K^2 wall
~K=48), masked block-Jacobi = today's solve, tau ridge in H and M, e/o Wiener statistic kept on the UNsupported x0 pair, solve_core
semantics, budget 2-3. G2 (coupled adjoint identity project_fplanes_mean_basis vs the coupled insert) is the one real risk.

See [[flex-pca-pcg-mstep-identity]], [[flex-pca-conventions-branch-state]], [[flex-pca-nuisance-deflation-flags-unset]],
[[flex-pca-solvent-mask-dominates-rna-basis]], [[flex-pca-must-be-knob-free]].
