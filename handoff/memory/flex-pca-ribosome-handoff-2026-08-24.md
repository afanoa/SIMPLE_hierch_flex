---
name: flex-pca-ribosome-handoff-2026-08-24
description: "Pointer to the EMPIAR-10028 ratcheting handoff doc: current best config, the DEFLATE_POSE patch, the spike benchmark, and the ranked next tests (COUPLED_DIAG is top)"
metadata:
  type: project
---

Full handoff: `doc/for_developers/ideas/flex_pca_ribosome_handoff_2026-08-24.md`
(companions: `flex_pca_10028_ratcheting_recovery_plan.md`,
`flex_pca_true_rotation_roadmap.md`).

One-line state: we DETECT the differential mode (latent z=+7.07, reproducible axis)
but basis capture is R2~0.16 against a benchmark-calibrated sortability threshold of
~0.3, so tails render only 0.5 deg instead of the motion.

**Top untested lever: `SIMPLE_COV_COUPLED_DIAG=1` was set in EVERY campaign run.**
It was adopted for a 1.68x speedup during the perf work and validated then on state
counts, never on basis fidelity; its own code comment says it drops cross-component
terms and is "off by default". The full per-voxel KxK coupled solve already exists
(`solve_coupled_basis_exp`). Test = unset the flag, score on the 2-deg spike and
real-data lobe R2. Cheap approximations if full is too slow: rotate-then-diagonal
(diagonalise the global latent covariance, solve diagonally in that frame),
orthonormal-basis EM, shell-wise KxK, one Jacobi sweep.

Best current config (reproduce in ~25 min on the 21k dev harness): bc64 lp16 neigs=20
n_probe_iters=25, `DEFLATE_BG=1 EM_DEFLATE=4`, **no** `DEFLATE_POSE`, cryoSPARC-pose
project `/usr/local/data/atifao/10028_dev20k/cs/cs.simple`.

Uncommitted code: `SIMPLE_COV_DEFLATE_POSE` in `simple_flex_pca_columns.f90` on
fwl-c142339 (backup `~/columns_backup_posedfl_applied.f90`, build rc=0). Keep the
code, keep the flag OFF -- it strips genuine differential signal because an off-centre
lobe rotation is not orthogonal to whole-map rotation. Patch text is in the handoff.

Benchmark harness (the durable asset): known lobe-rotation mode injected into real
images; ~10 min per estimator variant; calibration r=0.93@4deg, 0.75@2deg with good
basis vs 0.50/0.23 with POSE deflation on. Scripts: ~/spike.sh ~/inject.py
~/floor_eval.py ~/mstep_sweep.sh. See [[simple-star-pose-import-traps]],
[[never-signal-simple-processes-by-name]], [[flex-pca-must-be-knob-free]].
