---
name: flex-pca-ribosome-winning-recipe
description: "The config that produced visible 40S subunit rotation on EMPIAR-10028: full KxK coupled solve, BG+envelope deflation, NO pose deflation, 25 probe iters, cryoSPARC poses"
metadata:
  type: project
---

**Confirmed by eye on the delivered state maps (2026-08-24): the small subunit is
clearly rotating.** Full 105,247-particle run, merged 6 -> 5 states, mean-scale flat
(0.932-1.120), ~30 min wall at nparts=10 nthr=8.

## The recipe

Project: cryoSPARC-refinement poses, NOT our own alignment.
`/usr/local/data/atifao/10028_cspose/run_p1/cs/cs.simple` (105,247 ptcls; those poses
reconstruct to 3.42 A through our reconstructor vs 4.05 A for our abinitio3D).
vol1 + sigma2 from `1_bootstrap_rec3D/`: `recvol_state01.mrc`, `sigma2_it_1.star`
(`bootstrap_rec3D` generates sigma2 from the half-map difference -- plain
`reconstruct3D` dies on a missing `sigma2_it_0.star`).

```
simple_exec prg=flex_pca projfile=p.simple vol1=vol1.mrc \
    box_crop=64 box_rec=128 lp=16 mskdiam=320 ncols=256 neigs=20 \
    npreimages=6 n_probe_iters=25 column_sampling=lowfreq \
    state_axis=0 min_neff=2000 nparts=10 nthr=8 mkdir=yes
```

Environment:
```
SIMPLE_COV_EM=1                      # pure-EM path
SIMPLE_COV_EM_DEFLATE=4              # 4 envelope shells
SIMPLE_COV_DEFLATE_BG=1              # flat-in-mask background (solvent owned the 27x axis)
# SIMPLE_COV_DEFLATE_POSE            # MUST STAY UNSET -- strips the differential mode
# SIMPLE_COV_COUPLED_DIAG            # MUST STAY UNSET -> full per-voxel KxK coupled M-step
SIMPLE_COV_MERGE=1 SIMPLE_COV_MERGE_R=0.7 SIMPLE_COV_MERGE_LP=20
SIMPLE_COV_ONLINE=1 SIMPLE_COV_BATCH=40000 SIMPLE_COV_RELPRIOR=1
SIMPLE_COV_EM_MIX=16 SIMPLE_COV_POLAR_ESTEP=1 SIMPLE_COV_ZWHITEN=0
SIMPLE_COV_PROBE_CONV=999999 SIMPLE_COV_BASIS_MAX=200000 SIMPLE_COV_PROBE_MAX=200000
SIMPLE_COV_GPU=0 SIMPLE_COV_GPU_PREP=0 SIMPLE_COV_GPU_ESTEP=0
```

## The three settings that decide it

1. **`COUPLED_DIAG` UNSET** -- full K x K coupled solve (`rows=210 (full)`, 0.30 GB at
   K=20/bc64, no measurable slowdown). The diagonal approximation was in every earlier
   run because it was adopted for a 1.68x speedup during the perf campaign. User
   reports coupled is better across ALL their datasets; this ribosome run agrees.
2. **`DEFLATE_POSE` UNSET** -- measured: with it on, an injected 4-deg mode's basis
   capture collapses 0.46 -> 0.015 and sorting r 0.94 -> 0.50. An off-centre subunit
   rotation is NOT orthogonal to whole-map rotation, so the hard projection eats it.
3. **`n_probe_iters=25`** -- science parameter, never a speed knob. At 3 iterations the
   mode never separates from the noise bulk.

Output: 5 states, populations 81476 / 11723 / 6305 / 5600 / 143 (77% in one state --
the split is very uneven; compare the big state against the mid-sized ones, ignore the
143-particle state).

## Benchmark caveat -- do not trust the spike scores over your eyes

The 2-deg spike benchmark scored the coupled solve marginally WORSE (R2 0.263 vs 0.321,
r 0.742 vs 0.753) while the real data clearly shows the rotation. Both the spike scoring
and the earlier `lobe_rot` metric use the same crude density-weighted 2-means lobe
template (53/47 split -- not the real subunit), so they mis-rank real localized motion.
Use `~/local_motion.py` (template-free: patch cross-correlation -> displacement field ->
rigid fit to the DISCOVERED moving body, with the static set as control) and visual
inspection. Diagnostic inversion worth remembering: **LOW rigid-fit R2 with a genuinely
static remainder = REAL localized motion; HIGH R2 with a co-rotating "static" set =
pose jitter.** See [[flex-pca-10028-finds-real-mode-not-ratcheting]].

Also patched 2026-08-24 (source only, not yet built): dedicated
`COV_CALIB_MAX_PTCLS = 20000` / `SIMPLE_COV_CALIB_MAX` for `em_calibrate_noise_prior`.
That master-only pass was uncapped by default (`COV_BASIS_MAX_PTCLS = 0`) and read all
105k particles single-node before the first distributed round (~5 min at 1/10 capacity).
Backup: `~/columns_pre_calibcap.f90`.
