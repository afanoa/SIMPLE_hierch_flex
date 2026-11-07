---
name: pole-merge-line-inside-closed-pair-2026-10-22
description: "r1 pole merge: z1/z3 closed endpoints are distinct fixed points (Spearman 0.93-0.98 to each other, 0.99-1.00 across runs); a 0.95 merge line sits inside that range so P varies 3/3/4 between runs. Kept apart (P=4) the z3 seed becomes the closed extremum (0.93) and the z1 seed an intermediate: extrema reproduce at P=3 and P=4"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T08:09:14.096Z
---

Tested 2026-10-22 04:05, f692st node r1. Seeds were the argmax of four pole posteriors from `factorial/f692st_r1_polesets`, with greedy
pole-seeded stage A run for script / port1 / port2. Output `rnasep_klem/repro_r1/landscape_seeded4.txt`; scripts `repro_seeded.sbatch`,
`landscape_seeded4.sbatch`, `arm_readout_any.py`.

**Results** (readouts are closed-arm / open-arm):
- **z3-closed seed:** becomes the closed extremum, 0.93 / 0.52-0.55 in every run (K = 3 closed: 0.90 / 0.47).
- **z1-closed seed:** grows from 2.1-4.5k to 17-19k and becomes an intermediate (0.63 / 0.45), not a second closed state. It is
  within noise script-port2 (0.99) but 0.83-0.87 against port1, whose seeds had 2.1k rows.
- **Open A:** sharpens from 1.08-1.16 to 1.20-1.26.
- **K = 3 maps:** they lie inside the K = 4 landscape.

The port's integrated port3 "non-reproduced" fourth state (2.1k -> 17.8k) is this intermediate. It looked unmatched mostly by
construction, because the other runs seeded only 3 states.

**Why:** a merge line inside the distinct-pole range makes P differ between runs, and that breaks run-to-run reproducibility even though
each P on its own reproduces.

**How to apply:**
- Merge poles only at the replicate level (same fixed point, Spearman >= 0.99, calibrated by restarts), or always merge same-side
  endpoints. Never use a line in the middle of the observed between-pole range.
- Before calling a seeded state non-reproduced, check that the other runs seeded a counterpart.
- The user decides the rule.

Related: [[tree-leaves-not-reproducible-across-runs-2026-10-21]], [[landscape-must-be-reproducible-and-splits-data-derived]].
