---
name: robust-frame-loop-removes-heterogeneity-pose-bias
description: "2026-10-02 toy3: trimmed-LS frame correction of latent bins fixes POSES (extremes 7.25 -> 2.35 deg); CORRECTED 2026-10-04: its 'maps 0.75 -> 0.89' was a frame artefact (core-aligned maps 0.980 before and after), see pose-state-coupling-toy-verdict"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T06:11:31.538Z
---

Scripts on batch2 `~/toyJ/`: toyj_robust.py (`<toy> <fixed-EM run> <out oris> 7 recon <current oris>`),
robust_loop3.sbatch; data t3s0.02 (toy3, one region swinging +-15 deg, rigid core, SNR 0.02, 40k).

Mechanism found by debugging (align_debug*.py): refine3D_auto's bias is ONE coherent rotation per swing angle
(volume-frame error coherence 0.97 at the extremes, 6.7-7 deg). Plain least-squares alignment of the true-pose
extreme map to the centre map reproduces exactly that 6.7 deg (the consensus compromise), and Tukey IRLS started
there stays put (85% inliers). Trimmed LS (sum of the 50% smallest squared residuals, global +-10 deg rotation grid
on a voxel subset, then Powell) picks the largest matching part: 0.45 deg on true poses, 6.06 of 6.7 deg recovered.

Loop: fixed-pose EM latent (whole-molecule basis, NO regions) -> 7 equal-count bins along its leading direction ->
reconstruct each bin at the current poses -> trimmed-LS align to the centre bin -> R <- R A (R A^T is the wrong
direction). Round 1: pose error by |theta| 1.61/1.83/2.01/2.37 deg (refine3D 1.40/2.73/4.82/7.25, oracle 1.44);
region-1 maps at the extremes 0.895/0.893 (refine3D 0.75, true poses 0.988); converged after round 1.

**Correction (2026-10-04, [[pose-state-coupling-toy-verdict]]):**
- The "maps 0.75 -> 0.89" numbers (toyj_maps) were scored in the TRUE frame. A bias that is coherent per
  conformation only rigidly misframes each level's map.
- Scored after a rigid alignment on the static core (`~/toyPS/maps_al.py`), the extremes are 0.980 at the
  refine3D poses, 0.981 after the robust frame, and 0.989 at true poses.
- The pose-accuracy gain is real; the map gain is not.

Dead ends on the way (all on toy3): z-marginal (SOLVAR-type) joint search, conditional reference, mixture, latent
grid (needed hand regions, 4.91 deg), variance-weighted alignment (diverged), auto support from basis variance
(biased poses make static regions vary as much as the moving one). See [[joint-pose-heterogeneity-toy-verdict]].
