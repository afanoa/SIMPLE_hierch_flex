# flex3d on the raw fission halves (side session, 2026-10-31): hand-over

Full write-up: /Users/atifao/SIMPLE/doc/implementation_notes/grfp_tseries_flex3d_2026_10_31.md (untracked, not committed).
Cluster outputs: grfp_tseries_exp/pilot/flex3d (new folder only). New scripts in grfp_tseries_exp/scripts: fx3d_fit.py,
fx3d_fit_arm.py, fx3d_analyze.py, fx3d_diag.py, fx3d_view.py, run_fx3d_fit.sbatch, run_fx3d_arm.sbatch,
run_fx3d_analyze.sbatch. No existing script was edited. 121 core-h on csbdevel; no jobs left running.

Verdict: flex_pca (master 6ddc4a0d2, only lp changed: 1.5 / 5 / 16 A) adds nothing the occupancy model misses.
- lp 1.5: flex's own axis gate drops 55 of 60 axes; cross-half principal cos <= 0.47.
- Site changes follow the ramp fades weakly (best cross-half r 0.27-0.47) and add nothing vs the published segment
  maps (partial r <= 0.12; ramp 0.42-0.72, flex <= 0.28).
- Latents: 92-98% fast. Half-shared fast latents (p1, p3 at 1.5 A, up to 0.79) are view-predicted at 0.4-0.5 by
  same-orientation frames far apart in time: shared-pose misfit, not structure.
- The only half-consistent slow trend (p1 at 1.5 A) is the same dissolution the occupancy model already times.
- Off-site: nothing atom-scale and time-specific; no outward motion at fading atoms. NP87 at 16 A has a
  half-consistent change in the outer shell (5-9 A out, anti-correlated with the fades), most likely shared slow
  background.

Traps:
- Default lp 16 on a 57 A box is a 3-shell band.
- flex_pca_coordinates.txt holds deconvolved latents, rank one at 5-16 A; the MAP latents are in the raw block of
  flex_pca_embedding.bin.
- This build's embedding holds contrast at 1, but the leak is negligible here (R² <= 0.01).
