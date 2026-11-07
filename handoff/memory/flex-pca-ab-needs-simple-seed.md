---
name: flex-pca-ab-needs-simple-seed
description: flex_pca IS bit-reproducible at nthr=1 once SIMPLE_SEED is set; unseeded runs differ because calc_pspec bootstraps sigma2 from a random 25k-particle subset
metadata:
  node_type: memory
  type: project
  originSessionId: ab3f09e2-88fe-4a1c-aadd-cc6c1d565d5b
  modified: 2026-10-24T05:06:19.568Z
---

Measured 2026-09-30 on batch2 (toy project, 40.5k ptcls, box 96, 4.5k-particle pindfile subset, neigs=4,
nthr=1 nparts=1, default flex_pca incl. paired fit, merge, polish, embed, gridded state maps; ~2 min per arm).
- Unseeded: two runs of the SAME binary differed from the first data-dependent log line (mean self-scale
  2.9234E-01 vs 2.9232E-01) and in 47 of 60 output files (state-map corr 0.9999998, basis corr ~0.9998).
- Cause: `ensure_canonical_sigma_state` (simple_commanders_flex_pca.f90) runs calc_pspec when the project has no
  committed sigma2 state; `simple_calc_pspec_strategy.f90` caps the bootstrap at
  SIGMA2_BOOTSTRAP_MAX_PARTICLES=25000 and `ran_tabu%shuffle`s the E/O candidates, seeded from /dev/urandom by
  `seed_rnd` (called from parameters%new). The cap applies to ALL project particles, not the pindfile subset.
- With `export SIMPLE_SEED=<int>` (read by `simple_seed_base` in simple_rnd.f90) all 60 outputs were
  byte-identical across two runs of one binary, and across base-vs-patched builds (flex-GPU removal).

**Why:** the 2026-09-20 correction in [[flex-pca-distribution-changes-results]] blamed threaded reductions on this
same 40.5k toy; the random sigma2 bootstrap was also in play there, so that attribution is unproven (threaded
determinism under SIMPLE_SEED was NOT tested).

UPDATE 2026-10-24, THREADED RUNS ARE NOT BIT-REPRODUCIBLE even with SIMPLE_SEED. hetero_tree roots (nparts 4-8 x nthr 6,
seeded) gave different flex and fflex outputs run to run in the SAME layout: eigenvalues off by 2e-7 relative, z by up to
5e-4. Five IgM 100k runs over four builds fell on exactly TWO outputs. Every other root step (mean, cheap, recB0, site,
seeds, split) was byte-identical, and all part nodes were the same class (x6342). Likely cause: OpenMP float reductions,
which libgomp combines in thread-arrival order (flex_pca_deconv num/den, flex_pca_crossfsc rsum), or dynamic-schedule
accumulators; not yet pinned. Harness: stree_port/nproj_diag2.py. See [[hetero-tree-4p5m-scaling-2026-10-23]].
PINNED + FIXED (test only, 2026-10-24): schedule(dynamic) per-thread sums in the probe E-step, the embed loop and the basis EM
calibration, plus the reductions in calibrate_noise_scale and crossfsc_harvest_h. The fix (static schedules, ordered per-thread
sums) makes nthr=8 reruns byte-identical (toy, IgM 20k 1-part and 2-part), with no slowdown. It awaits the user's review:
.claude/patches/flex_determinism_on_056895a5c.patch. Benches: stree_port/flexdet/stree_flexdet{,3,4}.sbatch.

**How to apply:** set SIMPLE_SEED in every arm of any flex_pca A/B, and remember job scripts that unset all
SIMPLE_* env vars also unset it (re-export after the unset loop). For code-equivalence checks use nthr=1 +
SIMPLE_SEED + a pindfile subset and `cmp` every output: bit-identical is achievable. Harness:
`/mnt/beegfs/elmlund/afan/nogpu_ab/seeded/nogpu_ab_seeded.sbatch`. Related: [[flex-pca-exactness-verification]],
[[flex-pca-ab-testing-traps]].
