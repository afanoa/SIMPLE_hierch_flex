---
name: saturation-fission-site-test-verdict-2026-11-07
description: "VERDICT 2026-11-07 (Part 5 of igm-refs): knob-free per-site test = saturation of cross-half steps between 4 score-quantile groups with in-memory fission; fixes both qtoyflex GAPw failures with true groups and passes IgM 10 % a1, but loses 3-state arms split in two (p3: 3/5 vs GAPw 5/5): a complement to GAPw, not a replacement"
metadata:
  node_type: memory
  type: project
  originSessionId: 77c8a7fa-5cd0-4244-b81c-f308ec6cd639
  modified: 2026-11-08T00:03:20.203Z
---

Handoff `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/handoff/igm-refs/VERDICT.md` section Part 5; outputs in results/part5/.
Script `sat5.py` + `sat5.sbatch` (env FIX = qflex | x10fix | x10null, OUTD) in `math_reworks/latent/igm_disc/sat/`, on the latent
framework's cached planes (sparse_ppca/cache, read only; a guard refuses to rebuild a mismatched cache).

- Test: y1 = y + Z, y2 = y - Z per Fourier sample (planes P +- |CTF| Z; no stacks). Score t = y1 log-likelihood ratio between the
  input groups' site-masked maps. 4 groups = the f-quantile boundary + each side's median. y2 ML maps, cross-half squared steps
  q_A, q_mid, q_B. DISCRETE iff z(q_mid) >= 3 and one side has r + 2 SE < min(1, (2w)^2), r = q_side / q_mid, w = side fraction.
- Calls:
  - qtoyflex closed vs rest DISCRETE (r_A .36, via the rest side; the closed side is unresolved);
  - arm continuum at the s1 median not discrete (1.37 / 1.24);
  - 10 % fix a1 DISCRETE (.73 / .72 truth, .79 / .73 round-1);
  - null a1 not discrete (1.39 / 1.15 angle, 1.50 / 1.24 round-1);
  - random split: no call (z -0.6; z 23 without fission).
- Traps:
  - An unequal split needs the bound min(1, (2w)^2): a small side's narrow quantiles make any continuum look saturated.
  - The 2-state posterior prediction is overconfident (.22 predicted vs .73 for a true 2-state arm): never use it in the call.
  - Without fission z_mid measures the imprint.
- Costs:
  - Fission halves the score's information (10 % a1 d' 2.27 -> 1.55, r .50 -> .75), so the 5 % arms may give no call.
  - qtoyflex's round-1 closed split gives no call, because round 1 misses closed.

- Addendum (19:02, on Part 4b's p3 round-1 labels):
  - It finds 3 / 5 IgM main blocks: it LOSES both 3-state arms split in two (a2 1.10 / 1.37, a5 1.52 / .87; GAPw .22 / .19).
    The middle state straddles the boundary, so the steps are collinear and equal, as in a continuum.
  - 0 / 3 extras, 0 / 13 null blocks.
  - qtoyflex holo / closed: no call, because round 1's groups (31 %, 42 %) put the jump inside a side.
  - End-peaked control (IgM null a1, 20x keep-probability fall): dense side r 1.04 / 1.08 +- .4, not discrete.
  - So it complements GAPw; it does not replace it.

**Why:** the tree needs a site acceptance test without GAPw's calibrated .28. It must also handle curved continua and minorities
next to a continuum (RNase P's closed corner).
**How to apply:**
- Use sat5 (fis) only as a check where GAPw misreads (a curved continuum split in two); keep GAPw for IgM-like arms.
- Measure the 5 % arms before relying on it.
- Size jobs to ~7 min and 12 GB for an IgM root, 1 min and 3 GB for qtoyflex.
- Related: [[per-site-discreteness-test-gapw-verdict-2026-11-07]], [[rnasep-closed-state-must-be-retained]],
  [[data-fission-merge-test-verdict-2026-11-07]].
