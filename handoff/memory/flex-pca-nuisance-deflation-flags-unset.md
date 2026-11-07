---
name: flex-pca-nuisance-deflation-flags-unset
description: "flex_pca defaults leave solvent+envelope nuisance in the basis; on EMPIAR-10028 rank-5 deflation removes 85.7% of basis energy vs 9.3% at the rank-1 default"
metadata:
  type: project
---

`SIMPLE_COV_EM_DEFLATE` and `SIMPLE_COV_DEFLATE_BG` are **both unset by default**, so
`vdfl=1` and the mean-shaped deflation runs at **rank 1** — it removes only the rank-1
consensus-shaped term `(a_i - 1) * T_i P(R_i) mu` that exists because contrast is pinned
at 1 (`COV_UNIT_CONTRAST = .true.`, a deliberate choice: subtracting `a_i*T*mu` would also
delete conformational signal parallel to `T*mu`).

Measured on EMPIAR-10028 (P. falciparum 80S, 83626 ptcls, box_crop=64, lp=20, mskdiam=320):

| deflation | iter-1 energy removed | steady state |
|---|---|---|
| rank 1 (default) | 9.26 % | ~5-6 %, flat for 25 iters |
| rank 5 (`EM_DEFLATE=4` + `DEFLATE_BG=1`) | **85.73 %** | 10.72 % by iter 3, decaying |

**~86% of the delivered basis was nuisance.** The flat ~5-6% at rank 1 is not a failure to
converge — the nuisance is regenerated every iteration by construction, so deflation is a
per-iteration filter (the code's own 10076 figures are 34.6% -> 16.2% at convergence).
Rank 1 simply cannot reach the *shell-dependent* envelope family (B-factor spread, ice
gradient) or the flat-in-mask background mode.

Symptoms of running with the defaults, all measured on the 10028 baseline:
- `eig[0]/eig[39] = 27.30`, matching the code's note that the background mode "owns the
  leading axis (27x) whenever the mask is wide enough to admit solvent" (mskdiam=320 on a
  ~300 A ribosome).
- latent `var/noise` **flat at ~2.7** for components 2-39, only component 1 above it (14.5)
  — i.e. one nuisance axis plus a noise plateau, no conformational direction anywhere.
- eigenvalues not sorted: 14 ascending steps, tail rising 216 -> 786.
- rigid-rotation template R2 = 0.068 (k=5), 0.099 (k=10), 0.602 (all 40) — the rotation is
  only reconstructable by borrowing the anomalous tail.

The deflation defaults were "validated on 10076 only so far" — a *compositional* dataset.
See [[flex-pca-transport-rank-is-not-the-barrier]] and
[[flex-pca-solvent-mask-dominates-rna-basis]] (same failure, RNA, found via mskdiam).

## Attribution (2x2, measured on 10028)

| run | EM_DEFLATE | DEFLATE_BG | rank | comp1_snr | ratio | R2_k5 | signedness |
|---|---|---|---|---|---|---|---|
| baseline | 1 | off | 1 | 14.48 | 5.36 | 0.068 | 0.157 |
| bg1_dfl1 | 1 | **on** | 2 | 6.80 | 2.91 | 0.156 | 0.085 |
| bg0_dfl4 | 4 | off | 4 | **21.60** | **9.16** | 0.061 | **0.302** |
| dfl4bg | 4 | **on** | 5 | 6.68 | 2.84 | 0.208 | 0.130 |

**`SIMPLE_COV_DEFLATE_BG` does ~97% of the work. `SIMPLE_COV_EM_DEFLATE` alone is
HARMFUL** -- without the background template it strips envelope modes and leaves the
solvent mode *more* dominant than the untouched baseline (ratio 9.16 vs 5.36,
signedness 0.302). Never enable EM_DEFLATE without DEFLATE_BG.

`SIMPLE_COV_EM_DEFLATE=8` (rank 9) scored 0.201 vs rank 5's 0.208 -- saturated at 4,
and all 9 modes survived Gram-Schmidt, so the 1e-3 norm threshold does **not**
self-limit. A "request a generous rank and let it auto-trim" auto-rule will not work.

Tightening `mskdiam` is NOT a substitute: with BG on, R2_k5 falls 0.208 (320 A) ->
0.131 (300) -> 0.026 (280). A spherical mask cannot beat the background template --
even at 240 A the mask is still 24% solvent because the ribosome is not a sphere,
and by then 32% of the rotation signal is clipped. Deflation absorbs 1.000 of global
scale and **0.000** of rotation, so it costs nothing that a mask would preserve.
