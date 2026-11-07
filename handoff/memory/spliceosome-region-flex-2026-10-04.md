---
name: spliceosome-region-flex-2026-10-04
description: "10180 46k region-block flex at lp 12 (lp 8 box 320 was I/O-bound, 57 GB/part); auto B=2, forced B=4; SF3b-like region r3 axis q13 strongest; per-state refine3D_auto + bias-subtracted 4-state ladders launched"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T05:59:30.772Z
---

**Location:**
- Root `/mnt/beegfs/elmlund/afan/10180/autostates_c46k/`, consensus `10180/consensus46k_refine3D_auto` (5.91 Å, box 320 @ 1.699 Å, mskdiam 400).
- Scripts in `~/autostates/`: `as_stage1g.sbatch` (env H W TAG MSKDIAM BLIST LP NEIGS), `as_screen.sbatch` and `as_n4drv.sbatch` (cluster drivers with no 2 h limit), `as_staterefine.sbatch`.

**Band choice:**
- lp=8 with neigs=30 wrote 57 GB per worker part file (570 GB per iteration) and was I/O-bound for 90 min in iteration 1. It was aborted to `whole_lp8n30_aborted/`.
- lp=12 with neigs=20 takes 49 min for the whole fit; each n=3 ladder about 6 min.

**Regions:**
- Eigengaps B2 0.165, B4 0.104. Auto picks 2; B=4 has the better within/between ratio (1.44).
- Variability of the old 5 flex states (`sf3b_states_topvar.mrc`): B4 region 3 holds 52%, region 2 31%. B2 region 2 holds 84%.

**Screen (B4, `screenB4/`):**
- Null bias is small here: coarse SNR about 20, vs 300–500 on RNase P.
- Kept: q1 (r1); q7 (r2, z 17.7/8.7, var 254 vs null 122); q13 and q14 (r3, q13 z 14.4, var 282 vs 147); q19, q20, q21 (r4).
- Cross-region coupling r ≤ 0.10, including q7 × q13 0.087, a possible helicase+SF3b concerted motion.

**Launched 2026-10-04 ~02:00:**
- `staterefine_q13/s1..s4` + `s0` (random subset of equal size): refine3D_auto on the 4 equal-mass states along q13 (ZC=17). Plan: local resolution in region 3 vs the random subset vs the consensus.
- `screenB4/n4` + `n4_debiased` via `as_n4drv.sbatch` (SORT_LP=12).

Related: [[relion-multibody-benchmark-values]], [[flex-state-debias-model]].

**Model check (5NRL = EMD-3683, this dataset; `~/autostates/model_regions.py`):**
- Fit cc 0.655, SAME hand as the model.
- B4 regions = RELION's 4 bodies:

| region | RELION body | contents |
|---|---|---|
| r3 | SF3b | SF3b 0.90, SF3a 0.99, U2 snRNA 0.97 |
| r2 | helicase | Brr2 0.94 |
| r4 | foot | Snu114 1.00, U5 snRNA 0.94, Prp8 0.28 |
| r1 | core | U4/U6 proteins 0.98, Prp8 0.50, other tri-snRNP 0.65 |

- B2 (auto): r2 = SF3b + SF3a + U2 + Brr2 (the two mobile bodies together, matching RELION's concerted PC2); r1 = core + foot.
- Axes: q13 = SF3b, q7 = Brr2, q19/q20 = foot, q1 = core.
