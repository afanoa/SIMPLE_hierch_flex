---
name: toy2-hinge-swing-test-system
description: "2026-09-20: toy 2 = the visible-motion test system for multi-body flex (regions swung +-15 deg about hinges, ~10 A rms / 19 A tips; 27-state and 125-state sets under rnasep_regions/toy2, toy2c; generators toy2_make_vols.py/toy2_sim.sh). Settled: large swings ARE observable with enough PCs; rank-4 varimax global delivers each region's 30-deg motion at 0.99 like the 4-PC block fit; curvature needs 12 global / 4 per region; 24 global hurts; whole-molecule 8 states match exact references at 0.96-0.98"
metadata:
  type: project
---

Dirs: `/mnt/beegfs/elmlund/afan/rnasep_regions/{toy2vols,toy2,toy2cvols,toy2c}` (volumes vol_NNN.mrc, marginals
marg_r{r}_l{l}.mrc, masks reg96_*/env96/labels96_B3, gt_states.txt, gt_particles.txt; the *vols dirs also carry
gt_particles.txt for the scorers). Runs `toy2/T2_{B3_n4,B3_n8,G4v,G12v,G24v}`, `T2_G4v_ax{1,2,3}`, `T2_G12v_ax{1,2,3}`
(varimax attribution pc01->r3, pc02->r2, pc03->r1), `toy2c/T2c_*`. Launch: `TOYDIR=toy2 ~/launch_toy.sh <name> "<args>" "<env>"`.
Scorers: toy_score.py (attribution, localisation, theta and theta^2 R^2), toy_axis_score.py (corr + scale-free recovered
amplitude vs marginals; toydir = the *vols dir), toy_state_score.py (exact references). Numbers: note section 7.16.

**Why it matters:** the first toy's 2-5 % motions were invisible in a viewer and its whole-molecule states were
placement-limited; toy 2 shows the states themselves (whole-map corr with the mean 0.86-0.94) and gives the rank
ladder: 4 per region or 12 global capture the swing's arc, the delivery never needs the arc (backprojection),
24 global columns over-parameterise (4 of 24 reproducible, leaks 0.8).

**How to apply:** use toy2 (27 states) for any multi-body / rotation / placement change; toy2c (125 states) for
continuum behaviour; amplitude 15 deg is the calibrated visible-but-linear regime (25 deg gave pieces that no
longer overlap). Download for viewing: `~/Downloads/toy2_swings_20260920/`. Related:
[[flex-pca-varimax-rotation-recovers-region-axes]], [[flex-pca-state-placement-needs-motion-frame]].

**2026-10-14 addendum (toyL):** the three B3 masks (`rnasep_regions/fin_B3_{1,2,3}.mrc`) PARTITION the whole molecule
(region fractions 1.9/1.6/1.6 % of the box = the 5.2 % envelope), so swinging all three leaves no rigid core: at
+-30 deg a neigh refine3D from 10 deg off ends 15 deg (median) / 23-25 deg (extreme states) from the truth. For a
RNase-P-like toy keep region 1 (42 % of the mass) fixed (`AMPS=0:30:30` in `toyL/run_toyL_gen.sbatch`) so the arms
swing on a rigid core; the generator accepts amplitude 0 for a region.
