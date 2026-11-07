---
name: path-model-pose-marginalised-reweighting-2026-10-22
description: "VALIDATED occupancy route (swing30 + discrete simsites, 7/7 basin shapes): OT-morph path between two maps + per-particle likelihood with Laplace pose marginalisation + smooth reweighting (lam 0.01-0.001 chosen on truth, tau = pixel split ~1); swing30: TV 0.118, 2/4 kT barriers read 1.24/1.93, flat floor 0.25, 5 % state 3.8 %; real RNase P r1: two-sided p(s) with a 3.5 kT sparse region, not a few-state switch"
metadata:
  node_type: memory
  type: project
  originSessionId: f0f8499a-d0a6-450b-9f93-483f98f27006
  modified: 2026-10-22T04:15:52.233Z
---

2026-10-22 follow-up to [[arm-landscape-deconvolution-truth-tests-2026-10-21]] (user: "build a physical path model (morph between
deposited closed and open conformers), score every particle with per-state local pose search, validate on swing30 before RNase P").
Scripts batch2 routes/toy55: land_morph.py (+.sbatch), land_pscore.py (+.sbatch), land_reweight.py + land_reweight2.sbatch (positional args),
land_rwreal.sbatch, land_rwtable.py, land_twostate.py. Outputs rnasep_klem/landscape/{morph_*, ps_*.npz, rw_s30_*, rwreal_r1_*}.
- **Morph (unbalanced OT displacement interpolation, lp 8, blur 2 A, reach 40 A, no hand-drawn bodies).** swing30 t0 -> t30 reproduces
  the true rigid swing: corr 0.94 mean / 0.90 min in the moving region, against a cross-fade 0.81 / 0.64. The RNase P 18 -> 49 path moves
  33 % of conf18's mass, 86 % of it transported, mean displacement 12.6 A.
- **Scorer.** Per-particle quadratic expansion in 3 rotations (Fourier-gradient central slice) + 2 shifts; prior 2.5 deg / 1 px;
  lmarg = Laplace evidence. Pixel-split temperature 1.00 (synthetic), 0.86 (real). 108k x 31 maps takes ~30 min on 8 cores.
- **The method on swing30.** lmarg at tau_split with lam 0.01 (TV-chosen) / 0.001 (barrier-chosen): TV 0.118-0.124; dw2 1.24-1.31
  (true 1.92), dw4 1.93-2.21 (3.86); flat false floor 0.25-0.28; 5 % state 3.8-3.9 %. The morph path scores like the true path.
  - Without pose marginalisation (l0): dw2 0.23, dw4 0.42.
  - Tempering by the oracle tau* makes everything worse.
  - Held-out likelihood picks lam = 0 (spiky).
- **Few-state test** (land_twostate.py; free-position best 2 / 3 path states vs the full path, held out). Synthetic gains over the
  best 2: dw4 0.046, dw2 0.064, broad 0.064, flat 0.107 nats/ptcl. Real r1: 0.287 (over the best 3: 0.100 vs flat's 0.072).
- **Real r1 (92.5k, 18 -> 49).** Closed side s 0-0.3 holds 63 %, the middle 21 % with G rising to 3.8 kT at s ~ 0.8, and the open end
  s >= 0.93 holds 16 %. The two half stacks agree to 2 %. SIMPLE maps of the posterior-s windows are in
  view_rflex/landscape_f692st_r1_path.
- **Path sensitivity (r1 against 18 -> 28 instead of 18 -> 49).** Closed side 88 % vs 63 %, open end 4.6 % vs 16 %, sparsest point at s
  0.75-0.8 in both (5.7 vs 3.5 kT raw), few-state gain 0.19 / 0.06 vs 0.29 / 0.10. Robust: closed-side dominance, a continuous tail and the
  sparse region at ~3/4 of the way. Path-dependent: heights and end populations, so always name the path.
- **Basin rule (land_basins.py, lam 0.001, 0.6 kT, half check):** right on 7/7 truth shapes over two systems. The second system is
  discrete simsites: arm closed 8 % with nothing in between gives 2 basins, mass 0.080-0.086 (true 0.080), barrier read 6.4 kT. The
  latent rules cannot see this state.
- **Whole junk-free set (536k, 18 -> 49):**
  - 29 % closed side (s <= 0.3), 24 % middle, 46 % open end (33-40 % in the last bin: open states collapse onto the path end);
  - 2 reproduced basins, 0.54 / 0.46, barrier 1.9-3.1 kT;
  - few-state gain over the best 2 / 3 states: 0.229 / 0.076.
- **TRAP: ordering imprints kinematics.** Windows along a linear coordinate look lockstep (rank 1) by construction; windows along an OT
  morph look tip-first by construction. Neither is mechanism evidence.
**How to apply:** for occupancy/landscape questions use this route (morph between known end states, lmarg, tau_split, lam 0.001-0.01) and
quote barrier heights as lower bounds (35-50 % attenuation at RNase P SNR). Pass sbatch arguments with commas positionally, not via
--export ([[sbatch-export-splits-at-commas]]).
