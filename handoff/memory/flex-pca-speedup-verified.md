---
name: flex-pca-speedup-verified
description: The 2026-08-06 flex_pca speedup splits into a bit-identical free part and two env-gated levers that alone change the science
metadata: 
  node_type: memory
  type: project
  originSessionId: e31f34db-3e71-4d09-ab25-71f10f48fa4c
  modified: 2026-08-06T14:37:00.200Z
---

Verified on IgG (100k ptcls) on 2026-08-06. The speedup is **two separable things**, and
only the second can change a result:

1. **Free, provably inert (~1.6x).** Bounds hoist, fused SNR accumulate-and-clear,
   particle-parallel SNR, batched state insertion — plus the whole rest of that remote
   tree (auto-settings, bagging, merge restructure). Running the patched binary with
   `SIMPLE_COV_BASIS_STRIDE=1 SIMPLE_COV_GRAM_SP=0` reproduced the pre-patch run
   **bit-identically**: every state corr 1.000000, every PC 1.0000, eigenvalue max
   relative difference **0.000e+00**.
2. **The remaining ~1.35x** comes ENTIRELY from the capped initialiser subsample
   (`COV_BASIS_MAX_PTCLS=80000`) and the single-precision `ssyrk` Gram. These change the
   basis and are the only things to argue about.

**The pipeline is deterministic run-to-run** — two runs of the same binary were
bit-identical (eigenvalues 7e-11). So any A/B difference is real, never thread jitter.
This makes a repeat-run the cheapest and most informative control; run it first.

**How to judge "did the maps change": compare cross-arm state correlation against each
arm's OWN even/odd half-map correlation.** At `ncols=64 neigs=16` the worst cross-arm
state match (0.836) beat arm A's own even/odd agreement (0.822), i.e. the arms differ by
less than the data resolves. Individual eigenvector directions rotate freely (mean
subspace cosine 0.857) — expected given [[flex-pca-latent-is-997-percent-noise]] — but the
span and the delivered maps hold.

Starved configs mislead: at `ncols=24 neigs=10 n_probe_iters=1` the PRE-patch run was the
degenerate one (97.8% of particles in a single state); the patched run gave 5 populated
states. Don't A/B on a config that collapses.

All 8 `simple_test_flex_pca` unit tests pass, including `packed CG solve` (Sigma to 1e-6),
which is the numerical core the `ssyrk` change touches.

## Two stride defects found 2026-08-06 (one fixed, one open)

`cov_stage_subsample` divides the cap by nparts, which is right ONLY for a per-part list.
The probe passes one; the INITIALISER call site passes the MASTER's global list, so the cap
got divided twice and the stride was inflated by exactly nparts -- **stride 13, fitting the
basis on 7694 of 100000, at nparts=10**. Severity scales with nparts, so the handoff's own
nparts=4 Ribosembly runs were affected. FIXED by passing nparts only when
`flex_pca_is_worker()`; both branches then derive the same stride, which is what keeps the
SNR round's column choice matching the COLS round's accumulation.

STILL OPEN: **integer-stride quantization**. A cap of 80000 on 100000 wants stride 1.25 and
gets 2, keeping 50000 -- 37% of the allowed budget thrown away. Bresenham-style exact-count
selection (take i when `floor(i*nmax/n) > floor((i-1)*nmax/n)`) fixes it, stays deterministic
and halfset-safe. Also still open: selection is a stride over PROJECT ORDER, which correlates
with micrograph/defocus and ignores orientation coverage -- and angular coverage is exactly
what conditions a 3D covariance estimated from 2D projections.

Note `COV_PROBE_MAX_PTCLS=25000` is an ABSOLUTE count, not the "25%" it looks like on the
100k IgG set; on Ribosembly's 335k it is 7.5%.
