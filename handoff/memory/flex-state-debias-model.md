---
name: flex-state-debias-model
description: 2026-10-03 model-based selection-bias subtraction for full-band flex state ladders (~/autostates/as_debias.py) removes ~95% of null bias; S4 q2 amplitude 0.40 of the biased maps = cross-band estimate 0.41
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T03:38:40.237Z
---

**Model:** in state i of a ladder sorted along direction u, the bias is kappa · gamma(k) · dz_i · u(k).
- `dz_i` = state's mean coordinate minus the overall mean.
- `gamma(k)` = per-shell response, calibrated on the region's 3 null ladders (spectrum-matched random directions). It is zero above the 8 Å sort band.
- `kappa` = mean null coordinate variance / axis coordinate variance (noise share of the coordinate).
- The basis column is Fourier-padded from box_crop 72 to the state grid (160), at unit norm.
- Ladders live in `rnasep_autostates/<DS>/n4/{real_q*,null*_q*}_n4`; de-biased output in `n4_debiased/`.

**S4 validation (n=4), end-to-end SNR <20 Å:**

| ladder | before | after |
|---|---|---|
| null region 1, leave-one-out | 534–767 | 29–59 |
| null region 2, leave-one-out | 305–393 | 16–25 |
| real q2 (kappa 0.64) | 1104 | 172 |
| real q1 | 1262 | 146 |
| real q7 | 504 | 58 |
| real q8 | 402 | 33 (marginal) |
| real q10 | 292 | 10 (not real) |

- At 15–8 Å the nulls drop from ~5 to ~0.3–0.5.
- After/before amplitude for q2 is 0.40, matching the dilution-corrected cross-band estimate (0.41). The true change is about 40% of the old biased maps.

**Residual:** about 5% of the null bias energy remains, likely the other columns' correlated noise. A next fix would add a cross-column term using the other columns' per-bin dz.

Related: [[flex-band-sort-removes-selection-bias]], [[rnasep-thirds-selection-bias-null]].
