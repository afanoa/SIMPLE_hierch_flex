---
name: rnasep-klem-mask-clips-closed-arm
description: RNase P closed arm reaches 114 A from the box centre (90th pct 96 A) - any 100 A viewing/scoring mask clips the helix tip; klem's --mskdiam only masks the CALIBRATION consensus, its references are unmasked (mask reruns were no-ops); use 116 A masks for sharpening/scoring, reconstruct3D mskdiam 240
metadata:
  type: project
---

The user saw "the helix is never long enough" in the closed-state maps. Measured on the docked deposited closed map
(EMD-74209, conformer 18): arm voxels reach 114 A from the centre, 10 % beyond 95 A, 7 % beyond 100 A; the dense
molecule reaches 114 A (half-box 117 A at box 160 x 1.464 A).

What clipped it: the sharpened downloads (100 A soft mask), the identification/occupancy scoring masks (100 A on the
docked maps), reconstruct3D mskdiam 220 (1 % of the arm). What did NOT: klem.py's --mskdiam applies only to the
consensus volume used for the calibration and the noise model (line ~394, `VC = fftc(read_mrc(a.cons) * msk)`); the
E-step references are unmasked coarse-grid backprojections. Proof: the toy tree rerun with --mskdiam 232 reproduced
the 200 run to the digit (parent 0.61, children 0.64). I first told the user the references were blind to the tip;
that was wrong and was corrected.

**How to apply:** score/sharpen/occupancy masks at 116 A (sharpen_matched.py --rmask 116, score_maps / closed_arm_cells
masks), reconstruct3D / refine3D mskdiam 240 for the closed state; no need to re-run the real tree for the mask.
The arm itself: pure simulated 18 at true poses measures 0.85 by the metric, at realigned poses 0.81; real closed
subsets 0.78-0.81 (see [[rnasep-klem-real-data-setup]]). Related: [[klem-kernel-latent-em-2026-10-15]].

**2026-10-16 correction, the toy's closed tip is clipped at the source:** `simulate_particles` masks the input
volume (`build%vol%mask3D_soft(params%msk)`, line 106 of simple_commanders_sim.f90) and run_toy55_sim.sbatch passes
mskdiam=200, so every toy particle of conformer 18 was projected from a map whose arm beyond 100 A radius is gone.
The "pure simulated 18 at true poses = tip 0.64" ceiling was that mask, not a pipeline property; never compare the
real docked tip against the toy, and re-simulate with mskdiam >= 232 if the toy must contain the closed tip.
The real docked set's tip 0.75 vs the deposited 1.00 is a genuine deficit.
