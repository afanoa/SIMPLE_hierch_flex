---
name: empiar-10076-a200-inputs-are-broken
description: "10076 2D seeds: a200 is broken (unconverged, zero shifts); flexcls/10076/a2d_new is converged BUT on the protein-DARK stack (L17Combine_weight_local.mrc); the right seed is 10076_flex2D/a200n/e10076n.simple (converged, _neg = protein-white stack)"
metadata:
  node_type: memory
  type: project
  originSessionId: 03e8b343-4213-4208-890c-56346e7de4b2
  modified: 2026-10-01T22:30:00.000Z
---

2026-10-01: the user noticed the 10076 class averages "look bad". Measured: `10076_flex2D/a200`
(flex2D-branch snapshot f48eb527e, ncls=200, mskdiam=280, shared memory nthr=48) stopped at
iteration 8 with CONVERGED .NO., the registered cavg stack is cavgs_iter005 (box 160), and all
131,899 ptcl2D shifts are exactly 0.0 while e3 was assigned. Every 10076 benchmark against the
published labels (parents at chance 0.359) and every flex/dmap split montage used this project.
Pristine import: `10076_flex2D/import/e10076/e10076.simple` (132k, box 320, smpd 1.31, deftab).
Replacement run: `flexcls/10076/a2d_new/e10076.simple` CONVERGED (27 iterations, cavgs_iter027,
box 320); use it for every 10076 2D benchmark from now on (`flexcls/tier4n.sh <name> <args>`).
On it, 3-way vs the published labels: parents 0.357, flex 0.397 (41 state-pure subclasses),
diffusion maps 0.363 (4); band reproducibility flex 0.62/0.47 vs dmap 0.43/0.24 at 60-30/30-15 A. Related: [[empiar-10076-deposited-class-validation]] (same defect noted 2026-08-17 on
jclean10076), [[flex-cls-split-state-2026-09-30]].

**2026-10-03 correction:** `flexcls/10076/a2d_new/e10076.simple` points at the ORIGINAL protein-dark stack
`EMPIAR_downloads/10076/data/L17Combine_weight_local.mrc`; SIMPLE's mass/seg centring assumes white protein. The converged
abinitio2D on the inverted stack is `10076_flex2D/a200n/e10076n.simple` (`..._neg.mrc`, 27 iterations, seed purity 0.385
vs a2d_new 0.357). A flex2D x600 run seeded from a2d_new scored only 0.392/NMI 0.032 (10076_flex2D/seed_new_x600) vs 0.519
seeded from a200n earlier; diagnostics dg_R1_old400 / dg_R2_fix400 / dg_R3_fix600 (all a200n seed) separate the causes.
