---
name: arm-landscape-deconvolution-truth-tests-2026-10-21
description: "occupancy-landscape tests on swing30 (imposed p(t) by subsampling truth): naive landscapes flatten barriers 60-85 %; a 1-D per-particle coordinate cannot resolve a 2 kT double well at RNase P SNR under any calibration/smoothing, the full likelihood over true path maps can (1.7 kT); pose error sets the calibration (tau* 2.5); map-level unmixing gets basin masses but fills barriers"
metadata:
  node_type: memory
  type: project
  originSessionId: f0f8499a-d0a6-450b-9f93-483f98f27006
  modified: 2026-10-22T00:40:11.239Z
---

Task 2026-10-21 (user): quantitative arm dynamics / occupancy per position for RNase P, Dashti-2014 style but validated. Scripts (new, batch2
routes/toy55): land_truthmaps.py, land_e0(.sbatch) / land_e0b (oracle), land_deconv.py (smooth CV deconvolution module), land_e1 / land_e1b /
land_e1c (truth-free coordinate), land_joint.py, land_occ.py (closed-site occupancy units), land_real (real nodes), land_rec.sbatch (SIMPLE
window maps), land_confmix.py / land_leafmix.sbatch (unmixing onto deposited conformers). Outputs rnasep_klem/landscape/. Results note:
worktree intelligent-rhodes-abc9a6 doc/implementation_notes/arm_motion_landscape_2026_10_21.md.
- Shapes are imposed on simswing30's 108k swing particles by subsampling on the true t. No new simulation is needed. Shapes: flat, broad,
  dw2 / dw4 (quartic, minima 5 / 25 A), shoulder (5 % at t <= 3). The truth maps must be soft-masked like simulate_particles (mask3D_soft at
  mskdiam 200: cosine 62.6-74 px).
- ORACLE (31 true maps scored with lsamp_score2): own map first 11.8 % at true poses vs 8.0 % at node poses; posterior-mean corr 0.83 vs 0.72;
  temperature tau* 1.35 vs 2.46. Pose error, not noise, sets the likelihood calibration.
- Naive landscapes (argmax or point-coordinate histograms) flatten barriers: dw2 1.92 -> 0.78 (true poses) / 0.38 (node poses); dw4 3.86 ->
  1.3 / 0.6. The unpenalised NPMLE recovers dw2 at node poses (1.74 [1.28, 2.92]) but not dw4 (1.79), and it is spiky (TV 0.45-0.6).
- Truth-free per-particle coordinate (LS position on the end-to-end difference, cross-half, pixel sets A/B): |corr| 0.62-0.69 (the lsamp
  score axis 0.65-0.66). Split noise c 0.82 x Fisher against oracle residual 1.24-1.33 x Fisher, on every shape (view gain + pose error are
  invisible to a pixel split). A calibration from windows selected on the same lam is biased by the selected noise: use held-out pixels.
- Map-level unmixing (window maps NNLS onto the true basis): window mean t accurate to 0.1-1 A. Basin masses are right (low 0.499 vs 0.500,
  minor 0.097 vs 0.090 on dw2) and flat stays flat, but barriers are filled (dw2 0.6-0.7, dw4 1.2-1.5). NNLS pulls wells ~5 A inward and puts
  spurious 0.05-0.25 weight on the extreme basis maps (absorbing blur).
- Occupancy-unit route (land_occ, held-out windows, closed-site readout): the window readout is a compressed linear image of the truth
  (0.05 + 0.74 x true); the window-mean reliability gives c_reg ~ split. The smooth cross-half fit picks lam = 0 even under a one-SE
  rule: held-out likelihood cannot see spikes. Smoothed fits of dw2 are flat plateaus, so a 1-D coordinate (|corr| 0.67) carries too
  little information. c_win (window spread minus composition spread) came out negative. The joint map + particle fit (finite-difference
  gradients) collapsed its calibration to a constant.
- Map-level class shares depend on partition purity: windows of a noisy coordinate mix closed + open, and NNLS calls the mixture
  "intermediate" conformers.
**How to apply:** report basin populations (robust) with bands. Do not quote barrier heights from naive, map-level or 1-D-deconvolved
landscapes at this SNR. Quantitative barriers need per-particle likelihoods against a PHYSICAL path model (cryoER-style), better poses
(per-state local search), and a smoothness prior chosen on truth, not by held-out likelihood.
Related: [[landscape-sampling-stage-2026-10-21]], [[flex-pca-deconvolved-mixture-works]] (split replicate also missed ~40 % there).
