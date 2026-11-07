---
name: flex-pca-mean-scale-is-the-first-diagnostic
description: "flex_pca's per-shell mean amplitude scale must be flat near 1.0; when it isn't, vol1/the project is wrong and every state map collapses to consensus"
metadata: 
  node_type: memory
  type: project
  originSessionId: 0e59e2de-7b03-46a6-b41f-0c4b783e5d34
  modified: 2026-08-06T19:42:04.256Z
---

The `>>> FLEX_PCA per-shell mean scale s(sh) and s(sh)/s_global (D5)` block is the single most
diagnostic thing in a flex_pca log, and it is cheap to read. **The ratios must be flat near 1.0.**

Good IgG run (2026-08-05): `s = 3.6961E+01`, all 21 shell ratios in **0.968–1.150**.
Bad rerun (2026-08-06): `s = 3.2362E+00`, ratios **−1.530, 10.428, 7.686, 4.434, 1.072, then
0.037–0.285** — sign flip at sh=0 and a 280x span. Identical flex_pca settings otherwise (same
particles, ncols=128, neigs=20, min_neff=100, lp=15, d_tilde=177, probe on 100k, points=20).

**Why:** a non-flat ratio means `vol1`/the estimated mean disagrees with the particle data in
spectral shape, so the mean is subtracted with the wrong weight at every resolution. The residual
then carries a large spurious low-rank direction, which becomes z1.

**How to apply:** diff this block first, before suspecting any code change. The downstream signature
is unmistakable and also log-visible — reduced covariance max eigenvalue inflates (8.0e3 -> 1.2e4),
z1's eigval RISES while its rho FALLS (2.9e3/0.75 -> 4.8e3/0.53), tied bandwidth grows ~50%, mean
active states per particle jumps (2.77 -> 9.43), and one state absorbs most of the data
(21% -> 66% occupancy at 0.993 support fraction). At that point every state map is the consensus and
there is no conformational swing left to see. Causes to check: `vol1` post-processed/sharpened/
masked rather than the raw consensus, `vol1` from a different project than `projfile`, a `smpd`/box
mismatch after `reimport_particles`, or a missing `mskdiam` (its absence moved msk_crop 18 -> 30).

Full write-up: `doc/for_developers/ideas/flex_pca_swing_regression_audit_2026-08-06.md`.
Related: [[flex-pca-ab-testing-traps]], [[igg-test-project-setup]].

---

**REFINEMENT 2026-09-14: the per-shell scale is APPLIED, not just reported, and a non-flat ratio is
usually benign.** `simple_flex_pca_em_mean.f90:234` -- "3-point-smoothed, clamped radial scale
applied to the mean (FT state), then re-expand". So the printed `ratio=` column shows the CORRECTION
being made, not residual error. Do not read non-flat ratios as "vol1 is wrong and state maps will
collapse" (I did, twice, and over-corrected both times).

**Measured: pproc vs plain recvol as vol1 is IMMATERIAL.** On RNase P
(`reextract/24_refine3D_auto`, box 160, smpd 1.464), after per-Fourier-shell amplitude matching --
which is what the code applies -- `rec_final_state01_pproc.mrc` vs `recvol_state01.mrc` correlate
**+0.998 inside r<90 A** (the molecule), 0.959 whole box, and only 0.556 in solvent r>90 A where
pproc is hard-masked (59.8 % exact zeros vs 0 %). Solvent is outside the PCG support anyway
(0.038 of the box), so the choice does not matter. Its mean-scale ratios read 3.33/1.26/1.43/1.21
at shells 0-3 purely because the shell scale was absorbing the B-factor sharpening.
**TEST TRAP:** a real-space radial rescale does NOT emulate this -- it gave a misleading 0.769.
The code's correction is per FOURIER shell; match in Fourier space or the comparison is wrong.
