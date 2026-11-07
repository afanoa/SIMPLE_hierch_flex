---
name: flex-nuisance-block-code-map
description: "Where a FIXED nuisance block (micelle/background templates T, prior Psi) plugs into flex_pca's EM without changing ncomp shapes: Schur-complement marginalisation between the E-step former and the solver, residual subtraction at one CPU site; polar former needs T in its bank (or force Cartesian)"
metadata:
  type: reference
---

Mapped 2026-09-29 on flex-softblock (paths src/main/flex/, short names em/estep/iter/solve/mstep/embed/pmerge/lops).
- **E-step per particle** (fit_estep_solve_stats estep:272-338): statistics G=<PU_q,PU_r>, b=<PU_q,y>, c=<PU_q,Pmu>,
  e_mm, myv, contrast a=clamp(myv/e_mm,.1,5); solve probe_solve_ecm (solve:205-232): A=(a^2/s2)G+diag(prior),
  rhs=(a b - a^2 c)/s2; MCFA (COV_EM_MIX=16) from iteration 2 uses a full tied prior (probe_solve_mix solve:239-302).
- **Marginalise T exactly, keep everything ncomp-shaped** (insert before the solver, estep:~283):
  S = G_cc + (s2/a^2) Psi^-1; G' = G_zz - G_zc S^-1 G_cz; b' = b_z - G_zc S^-1 b_c; c' = c_z - G_zc S^-1 c_c;
  E[c|y] = S^-1(b_c/a - c_c - G_cz E[z|y]). Needs per-particle P T_j (projections of each template).
- **Formers:** polar is always on (l_pol_es=.true. iter:1207; bank fit_polar_bank_build estep:32-184 holds [mu;U] as
  slots 0..ncomp) -> T must be appended as bank slots, or force the Cartesian former (on-the-fly
  project_fplanes_mean_basis lops:463-611; much slower).
- **M-step residual:** one CPU site, in place, estep:326-335 (subtract_mean_banded / full-plane) -> subtract
  a*sum_j E[c_ij] P T_j there (plug-in; drops the Cov(c,z) term). GPU/bank residual sites iter:735-738, 822-824, 915
  must be gated off. Prior-variance update from posterior second moments (estep:313-315 -> mstep:58-60, :430).
- **T must stay OUT of ncomp/basis_recs** (orthonormalize_blocks basis:110-168, flex_col_block pcg:2678-2691, merge).
  Load T+Psi at both fits iter:1462-1488, worker iter:1714-1721, polish model:344, merged embed embed:168-181.
  Reuse compose's load/pad loop (compose:86-144) under a separate namespace; do not call compose_basis_from_runs
  (clobbers namespaces, throws <2 cols, diverts from EM). Keep U orthogonal to T via deflation templates
  (mstep:272-390, merge_deflate_mean_shaped pmerge:971+).
Limits: neigs <= 48/region; device E-step ncomp+1 <= 64. Side bugs flagged as a task: merge ignores DEFLATE_BG/
DILATION opt-outs (pmerge:986,993); l_gpu_probe lacks the pcg guard (iter:546). Related: [[micelle-marginalization-proposal]].
