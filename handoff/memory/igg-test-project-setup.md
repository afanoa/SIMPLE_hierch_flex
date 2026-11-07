---
name: igg-test-project-setup
description: How to make the IgG heterogeneity dataset runnable for flex_covariance testing
metadata: 
  node_type: memory
  type: reference
  originSessionId: 5aaa37ef-99ec-4ef0-95d1-add578e1f48f
  modified: 2026-08-03T15:55:50.991Z
---

IgG simulated heterogeneity dataset (SNR 0.01, 100k particles, box 128, 3.0 A/px,
50k/50k even/odd). Refined poses + consensus map live in
`~/manifoldem/remote_igg/IgG-RL/igg/4_reconstruct3D/igg.simple` (consensus =
`recvol_state01.mrc`).

The project's stk paths point to a remote prefix `/usr/local/data/atifao/IgG-RL/...`
that does NOT exist locally (and /usr/local/data is not writable). Fix by repointing
to the single merged local stack (same particle order, poses preserved):

```
cp 4_reconstruct3D/igg.simple <workdir>/
simple_exec prg=reimport_particles projfile=igg.simple \
  stk=~/manifoldem/remote_igg/IgG-RL/images/snr0.01/igg_manifoldem_merged.mrcs ctf=yes
# -> writes <workdir>/1_reimport_particles/igg.simple with local paths + refined poses
```

The merged stack `igg_manifoldem_merged.mrcs` is 128x128x100000 float. Do NOT
re-import from the star files (snr0.01_simple_paths.star / igg_manifoldem_merged.star)
— those carry the ORIGINAL unrefined MD-frame poses, not the SIMPLE-refined ones.

Run: `simple_exec prg=flex_covariance projfile=1_reimport_particles/igg.simple
vol1=<consensus.mrc> box_crop=64 lp=15 ncols=24 neigs=6 npreimages=5 nthr=8`.
Use fromp/top to subset. See [[flex-covariance-recovar]].

Fresh binary after builds is `build/production/simple_exec` (build/bin/simple_exec
can be stale). No sigma2 star present -> flex_covariance falls back to unit
white-noise spectrum (correct for this synthetic white-noise data).

**Per-particle latent scoring cannot be run on this laptop** (confirmed 2026-08-03).
`scripts/score_latents.py` — the instrument for the metrics that still matter
(`cv_R2` of GT-PC1 from per-particle z, `pp_best`, `frame_capture`) — hardcodes
`/home/afan/simple_igg/analysis/gt_frame_coords.npz`, which is not present locally
and has no local equivalent. Only `flex_basis_capture.py` runs here, and that metric
is RETIRED (see [[flex-covariance-basis-defect]]), so a local run can confirm a
config executes and time it, but cannot rank configurations. Do estimator
comparisons on the machine that has the GT frame coordinates.

Run outputs land directly in the CWD with `mkdir=no`; the program and its file
prefix were renamed `flex_covariance` -> `flex_pca` (`flex_pca_pc*.mrc`,
`flex_pca_state_*.mrc`, `flex_pca_coordinates.txt`). `flex_pca_mean.mrc` is NOT
written in the default vol1 path (`COV_MEAN_FROM_DATA=.false.`), which is why
flex_basis_capture.py's frame control reports `nan`.
