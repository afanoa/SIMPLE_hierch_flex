---
name: flex-pcg-cold-start-and-fsc-prior
description: "2026-09-12: the flex PCG M-step now starts from ZERO (b12d8f6bd, budget-0 parity 1.6e-7); the FSC prior was ALREADY in the operator (set_ridge) but is inert in every single-fit run because only the paired engine writes the cross-FSC records; when active it is 0.52 of the data term and still does NOT make maxits_pcg=2 sufficient (residual 0.81 vs 0.78 without it), so production's iteration count does not transfer to flex"
metadata: 
  node_type: memory
  type: project
  originSessionId: e5141833-c145-408a-ba48-7c3a3ffe0fbd
  modified: 2026-09-12T18:28:57.555Z
---

**Cold start (b12d8f6bd, branch flex-pca-pcg, batch2 snap_be8c70414).** Every positive budget now runs CG
from zero, the production reconstruct3D_pcg convention; budget 0 stays the explicit gridding mode and ships
`solve_coupled_basis_exp` untouched. Parity gate: budget-0 eigenvolumes vs run 56 differ by 1.6e-7 (float32
storage noise), eigenvalues bit-identical, components 5-8 exactly equal. Rationale for dropping the warm
start: `solve_coupled_basis_exp` IS this operator's block-Jacobi preconditioner applied to the numerators, so
a cold CG's first iterate already lands there; the warm start bought an off-by-one, cost the A-norm
convergence bound (stated relative to x0) and carried an unprojected object into a projected solve.

**The FSC prior was already implemented in the operator** -- I wrongly claimed otherwise before checking.
`apply_operator` carries a per-component, per-shell Fourier ridge (`self%ridge`, `l_ridge`), installed by
`set_ridge(fit%xf_invtau2)` at simple_flex_pca_em_mstep.f90:100, and the same invtau2 goes into rho via
`add_invtausq2rho_coupled` so the preconditioner and the gridding solve agree with the operator by
construction. `xfsc_build_invtau2` derives it from the cross-FSC with a lag (`crossfsc_latest_upto(it-1)`),
two arms (1 cross, 2 = per-shell min(internal, cross)).

**It is inert in every single-fit run.** `ctx%l_writer = l_paired .and. l_master`: only the PAIRED engine
writes the cross-FSC records the ridge consumes. A single fit with `SIMPLE_COV_XFSC_REG=2` logs
"XFSC ridge requested but no artifact exists and no writer is active; it stays inert" and
"arm degraded to 0: no record stamped <= t-1". So NO PfCRT result any of us has looked at had an FSC prior
anywhere. Paired needs `SIMPLE_COV_PAIRED=1` plus either `SIMPLE_COV_PAIRED_MERGE=1` or
`SIMPLE_COV_PROBE_ONLY=1`, and refuses the complement block.

**Scale audit (be8c70414, `SIMPLE_COV_PCG_VERBOSE=1`)** reports each regulariser against the data term on the
rhs probe. Measured: Tikhonov/data = 1.2e-3, matching `lam_rel=1e-3` exactly, so lam is correctly scaled
against the doubled-lattice Toeplitz operator (no padf^3 error, which I had suspected). With the prior active
from iteration 2: **FSC-prior/data = 0.52**, and lam doubles because it is derived from rho which now carries
invtau2.

**The budget does NOT transfer, and the prior does not rescue it.** Cold residual after 2 iterations on PfCRT
(box 100, protein support, 1.7% occupancy): 0.78 without the prior, 0.81 with it -- unchanged. 20 iterations
reach 2.8% (8.8 s/solve vs 1.7 s; probe iteration 70 s vs 29 s). Two iterations can only build a 2-dim Krylov
space for a coupled KxK system on a 2%-occupancy support; production gets away with 2 because it solves one
volume on a full sphere. **Conclusion: convention parity is achievable (cold start, prior in the operator) but
`maxits_pcg=2` is not scientifically justified for the flex M-step; flex wants a stopping rule (rtol>0) with a
cap, not production's number.** Runs: 91 (budget-0 parity), 94/95/96 (cold ladder 0/2/20), 103/104
(prior requested, inert), 107 (paired + prior), 108 (paired control).

**Defaults aligned (7c327db3b, snap_7c327db3b):** the COMMANDER now matches the UI declaration,
FLEX_PCG_MAXITS_DEFAULT=20 and FLEX_PCG_RTOL_DEFAULT=1e-3, so a UI-driven launch and a bare command line run
the same estimator (finding F7 closed). rtol>0 matters beyond the residual test: it ARMS the
FLEX_PCG_XTOL=1.5e-2 diminishing-returns stop on dx/x (`stop_xtol = rtol > 0.0 .and. dxx <= XTOL`), so 20 is
a cap with two stop rules rather than a fixed count. On PfCRT neither fires (dx/x 3-4% at 20, stop=maxits).
Verified by running with no keys on the command line: iters=20, stop=maxits, |x0|=0. NOTE `script=yes` writes
the USER's command line, not the derived defaults, so the script is not a way to check them.

**State reconstructions** take `max(maxits_pcg, FLEX_PCG_MAXITS_FLOOR=5)` cold on the SPHERE, so they now get
20 as well. Measured from flex_pca_state_pcg.txt (the per-(state,half) table every states run writes):
5 iterations leave 0.069-0.092 residual / 6.9-11.6% update; 20 leave 0.033-0.036 / 3.1-4.0% for populated
states (~4000 ptcls) and 0.060-0.067 / 8% for near-empty ones (38-46 ptcls, prior-dominated anyway). The
floor of 5 now binds only when maxits_pcg=0 is passed to select the gridding basis. Box-300 state solves are
the dominant cost of the state stage and scale with the count.

See [[flex-pca-pcg-branch-state]], [[flex-pca-complement-block-and-step0]], [[pfcrt-injection-floor]].
