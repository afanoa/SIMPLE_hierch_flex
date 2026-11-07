---
name: motion-series-prototype-2026-10-24
description: "Phase-2 motion series (5 posterior bins per site -> one 5-state reconstruct3D on the tree's root stage-B project): smooth monotone series, IgM ends ~9.6 A at 10k; OLE at box_crop 200 ~5.5 A per 16k bin; recipe, costs (5.1 core-h) and the 2-column posterior TRAP"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-24T22:56:34.743Z
---

2026-10-24, the coordinator's phase 2 of the user's "extremes + how we go between them" task; tools in routes/toy55
(motion_series2.sbatch, motion_bins2.py, motion_score.py, simple_fsc.py).
- Recipe, no new code:
  1. cp node_r/stageB/b.simple; update_project qsys_name=local.
  2. `selection oritype=ptcl3D infile=states.txt` (0 = out, 1..5 = bin).
  3. The tree's own reconstruct3D line from node_r/recB0/hetero_step.sh, with nstates=5 and postprocess=no (unsharpened, which skips
     the costly per-state postprocess) [+ box_crop].
  4. SIMPLE_QSYS=local on snapshot m8c52t.
- Costs: an IgM 200k box-160 rec3D takes ~110 s on 12 cores. OLE 88k at box_crop 200 takes ~130 s on 16 cores. OLE at its native
  box 360 would take ~30 core-h per 5 states (assembly + postprocess dominate).
- Results (genuine-quantile bins; junk kept in the same score intervals):
  - IgM M r_s1 / r2_s4 and m8c17 r2_s1: truth means even and monotone; neighbour corr 0.99+; ends 0.93-0.97; FSC 7.1-7.6 A at 32k
    per bin, ~9.6 A for 10k end bins. M's middle bins are ~45 % junk but still resolve, because contrast weighting down-weights it.
  - OLE r1_s4 / r2_s4 / r_s2: 5.1-5.7 A per 16.1k bin; ends corr 0.87 / 0.88 / 0.78. Maps are in ole_12707/motion_series/<site>/
    (bin1..5 + halves/); delivery is the coordinator's call.
- Promoting the weak IgM motions' deep candidates to the root gives no lift (span ratio 0.43-0.51, abs(rho) 0.44-0.52). An oracle
  removal of truth junk adds only +0.01-0.02, so junk is not the limit.
- Residualising the weak score on the other population-wide kept-site logits (cond_resid.py) lifts only a3_th (0.45 -> 0.58,
  ratio 0.44 -> 0.57); a5_th and a5_ph do not move or get worse. Regressing out root fflex z also removes the target (a3_th -> 0.21);
  don't include it by default.
- Conditional in-branch series are clean (0.62-0.70 of the level oracle) but small (4-8k). For each weak motion, 3-5 disjoint branches at
  abs(rho) >= 0.5 cover 68-88 % of the population; branch-level scores of 0.51-0.58 exist over 65k-95k particles.
- TRAP: factorial post files can hold one column per pole (m8c17: 2 columns). Reading them flat doubles the rows; motion_bins2.py takes
  the pole-1 column.
**Why:** this is the deliverable design for motion extremes + paths: one continuous score per motion, binned and reconstructed.
**How to apply:** reuse motion_series2.sbatch; price OLE-sized data with box_crop and postprocess=no; check each post file's column
count. Related: [[igm-motion-extremes-capture-2026-10-24]], [[discrete-splits-capture-the-span-states-that-matter]].
