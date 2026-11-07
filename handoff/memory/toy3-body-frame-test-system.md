---
name: toy3-body-frame-test-system
description: "toy3 = one RNase P region swinging continuously (31 levels, +-15 deg) with real fine detail; the 10-minute loop that tests per-body poses (body-frame reconstruction) against binning; results 2026-09-29"
metadata:
  node_type: memory
  type: project
  originSessionId: a0f08b9d-781d-471c-b0a6-698a8a8afe30
  modified: 2026-09-29T18:51:46.021Z
---

Built 2026-09-29 as the test "the current way fails, body-frame poses might work". Dirs on batch2:
`/mnt/beegfs/elmlund/afan/rnasep_regions/toy3/{vols,s0.05,s0.02,logs}`, scripts `~/toy3/` (also in the session
scratchpad): `toy3_make_vols.py` (region 1 of fin_B3 swung about a hinge at its inner boundary, 31 levels over +-15 deg,
regions 2/3 static; writes piece_r1_l0.mrc = exact unrotated body, bodymask_r1.mrc, static.mrc, hinges.txt in centred xyz
voxels), `toy3_body.py` (composes per-particle poses with a rigid body motion: M = R_i Q, d = R_i(c + t - Q c); modes gt /
gt cen / latent / bin), `toy3_score.py` (masked FSC vs exact references; BIN/FLEX modes build the weighted reference from
gt_particles; `sharp` = nearest true volume), `common.sh` (`rec <dir> <oritab>` = new_project + import_particles +
optional `selection infile=<tab>.flags` + reconstruct3D, ~40 s for 40k box-96 particles), `arm.sbatch` + `loop.sh <SNR>
<TAG> <MODE> "<ARGS>"` = one arm per 16-cpu job (compose, rec, score -> s<SNR>/SCORES.txt, ~90 s). Block fit
(`neigs=2`, labels96_B3, nparts=1 nthr=16) takes 3 min; state_axis ladder resumes 30 s each.

**Conventions found:** SIMPLE image = P_z[V(R^T x)] with R = euler2m (row-vector matmul), so a body rotation Q composes
as M = R_i Q. The gt rotation needs osign=-1 relative to the generator's Rodrigues in zyx index space; the latent
generator projection needs osign=+1 (its sign is opposite). In-plane shift sign RESOLVED 2026-09-29 on `toy3_offc` (HINGE_MODE=outer, hinge at the region's outer
boundary, |d| up to several px): ssign=-1, i.e. new x,y = old x,y MINUS the projected body displacement d (gt arm
o-1_s-1 = 0.988 at 8 A; s+1 0.36, s0 0.64). The base toy's hinge at the box centre made d ~ 0. `import_particles` forces state=1 (line 274 of
simple_commanders_project_ptcl.f90); exclude particles with `selection oritype=ptcl3D infile=<0/1 per particle>`.
numpy 1.14 on batch2: no np.quantile (use percentile), no fromfile(offset=).

**Results (SNR 0.05, 40,300 particles, true poses; metric = masked FSC in the body at 8 A / 6 A vs the exact unrotated
piece):** consensus floor 0.63/0.57; body-frame with TRUE angles 0.97/0.95 (Nyquist 4.88 A); latent-derived motion (6
rigid generators fitted to U z_i in the body mask) raw 0.71/0.63, unsupervised auto scale (match U(zbar_N - zbar_1) to
the on-axis end states, scale 2.9) 0.89/0.81, oracle scale 4.45 0.94/0.88. Latent on-axis rotation vs true angle: corr
0.96, amplitude compressed 4.45x (MAP shrinkage + linear regime), off-axis rms 0.7 deg. SNR 0.02 is indistinguishable at
full N (floor 0.63, auto 0.88). Core cost: the static regions in the body-frame map drop from 0.996 to 0.78-0.89 at 8 A
-> composite (core from the consensus, body from the body-frame map). On-axis ladder states reproduce their own blurred
references at 0.95-0.99 even at 24 states / 2.7k neff, so "vs exact reference" does not measure sharpness; use the
sharp-target metric.

**How to apply:** iterate the latent->motion mapping with `loop.sh` arms (90 s each); the open items are the amplitude
calibration (auto recovers ~2/3), the shift convention (off-centre hinge toy), and the RNA transfer. Related:
[[toy2-hinge-swing-test-system]], [[multibody-flex-means-region-pca]], [[rnasep-combined-block-embedding-states]].

**Low-N variant (`toy3_lowN/s0.01`, 12,400 particles, SNR 0.01) = the test the current way fails.** FSC of the moving
domain vs the SHARP true volume at 8/6 A: consensus 0.59/0.51; oracle 12 bins 0.65-0.76/0.40-0.61; flex on-axis 12
states 0.73/0.50; best flex slice (6 states) 0.83/0.70; body-frame latent uncalibrated 0.65/0.56; global scale 4-5
(oracle) 0.88/0.77; TRUE angles 0.94/0.88. Base toy (40k, SNR 0.05) is too easy: slicing ties body poses.

**Amplitude calibration, resolved 2026-09-29 evening.** The 4.5x compression is NOT tangent saturation (25 % at 15 deg,
`curve` mode), NOT band mismatch (band-limited generators only 4.5 -> 3.4, LATENT_LP env) and NOT MAP shrinkage (cached
posterior precisions give 1.15x, `deshrink` modes). What works, label-free: `autog:<ladder 1_flex_pca>` = project BOTH the
on-axis end-state difference (state_N - state_1, image averages, amplitude-faithful, greyscale-matched to vol1) and the
latent difference U(zbar_N - zbar_1) on the six rigid generators and take the ratio of coefficient vectors: scale 4.79
(cos 0.999), slope to truth 0.97, FSC 0.885/0.771 = the oracle optimum. The direct map-to-map version (`auto`) gave only
2.9 (cc 0.69) because the state difference carries non-rigid content. Sharpness by band power is NOT a usable
label-free criterion (monotone in scale). Recipe per region: block fit -> 3-state on-axis delivery -> `autog` -> compose
-> reconstruct3D -> composite with the consensus core.

**2-DOF variant (`toy3_2dof/s0.05`, region 1 swings about TWO perpendicular axes, 11x11 levels +-15 deg, 40,414
particles, SNR 0.05; generator env `AXIS2=15,11`, the second angle rides in the theta2 slot of gt_states) = THE TEST THE
CURRENT WAY FAILS.** vs SHARP truth at 8/6 A: consensus 0.43/0.37; oracle bins along axis 1 (6 or 12) 0.59-0.62/0.51-0.59;
flex on-axis ladder 3/6/12/24 states 0.48-0.65/0.44-0.58 (every 1-D slice leaves the other axis blurred); body-frame
latent (map-to-map `auto`, scale 2.8) 0.84/0.74, uncalibrated 0.53/0.42; TRUE angles 0.983/0.973 (both rotation
composition orders identical). Latent recovers both axes (corr 0.92/0.93 with the two true angles) but with a spurious
third-axis component (off-axis rms 6 deg at the calibrated scale): restricting the rigid fit to the dominant rank of the
omega covariance is the next lever. `autog` arms on this variant submitted 2026-09-29 ~19:00.
**2-DOF autog result:** scale 3.75-3.84 (cos 0.999), slopes to truth 0.91-1.01 on BOTH axes, FSC 0.858/0.762 (ceiling
0.983/0.973). The "off-axis rms 8.7 deg" diagnostic is misleading in 2-DOF (it measures everything not on axis 1); the
fitted rotation covariance is exactly rank 2 (singular values 2.49/2.32/0.00 deg), so `LATENT_RANK` restriction changes
nothing. The residual gap to the ceiling is per-particle noise in the two angles (corr 0.92 each), i.e. the next lever
is a per-body local pose refinement (tier 2) or a denoised latent, not the mapping.

**Ribosome 10028 transfer (2026-09-30, `10028_simple/RIBO`, scripts `~/toy3/ribo_{setup,ladder,arm}.sbatch`):**
co-variation split of the 20-column W_refit basis is weak (ratio 1.2; pieces 49/51 % of the envelope, piece 1 holds
0.55-0.71 of every column's energy and is the soft half: 3.4/4.0 vs 3.3/3.6 A). Calibrated motions 1.1 deg (piece 1)
and 0.3 deg (piece 2) median, autog cos 0.94 / 0.72. Against the RE-GAUGED control (see
[[reconstruct3d-negative-e1-degrades-map]]; the raw control at 4.5 A is an artefact): piece-1 body frame 3.6/4.3 in its
mask vs 3.4/4.0 (worse), piece 2 unchanged. NEGATIVE on real data, as on RNase P. Both real datasets: the latent carries
no rigid body motion large/clean enough for per-body poses to beat the consensus; the toy line remains the only
positive evidence.
