---
name: flex-pca-per-shell-snr-unconfines-the-basis
description: "Normalising the flex_pca SNR proxy per shell moves the basis from 115 A to 44 A and unconfines it from >60 A, but halves cross-halfset reproducibility — signal-vs-noise not yet settled"
metadata: 
  node_type: memory
  type: project
  originSessionId: 991559ff-d34f-451c-ad5a-ac1ba95b88c7
  modified: 2026-08-10T17:14:50.618Z
---

`SIMPLE_COV_SNR_SHELLNORM=1` (added 2026-08-10, default OFF) replaces the SNR proxy's single global
noise floor with a per-shell one: `snr = max(0, (var/dens)/q25(sh) - 1)`, i.e. excess over the
shell's own level, in units of it, so voxels at different resolutions are finally comparable in the
greedy argmax. q25 rather than median because the median of a signal-rich low shell is a signal
estimate, and normalising signal by signal would penalise exactly the shells the conformational
covariance lives in.

Measured, IgG 30k, box_crop=64 lp=12, nparts=5, bagged:

| | >60 A | 60-40 | 40-30 | 30-24 | basis res | mean cos | top1 |
|---|---|---|---|---|---|---|---|
| default (global floor) | 93.11% | 1.99% | 4.08% | 0.32% | 114.7 A | 0.5192 | 0.868 |
| per-shell | 27.71% | 19.54% | 35.82% | 16.79% | 44.4 A | 0.2431 | 0.607 |

Columns move out from |xi| mean 8.0 -> 12.5. So it DOES what it is meant to — the basis is no longer
a >60 A object — but cross-halfset reproducibility halves.

**Do not conclude "the high-resolution content is noise" from that alone.** Reproducibility is
exactly the metric a smoothness bias would corrupt, and a 44 A basis vs a 115 A basis is the worst
case for it. The one point against pure smoothness bias is that the same metric detected real
information loss at lp=36 (see [[flex-pca-basis-lives-above-60-angstrom]]).

**THE SCALE-N TEST SAYS IT IS SIGNAL.** Real signal reproduces better with more data; noise does not.
Run at 100k, same day:

| config | N | mean cos | delta with N | basis res | >60 A | 40-30 A | 30-24 A |
|---|---|---|---|---|---|---|---|
| default | 30k | 0.5192 | — | 114.7 A | 93.11% | 4.08% | 0.32% |
| default | 100k | 0.7010 | +0.182 | 148.6 A | 99.11% | 0.36% | 0.05% |
| per-shell | 30k | 0.2431 | — | 44.4 A | 27.71% | 35.82% | 16.79% |
| per-shell | 100k | 0.4799 | **+0.237** | 61.3 A | 58.54% | 28.81% | 12.11% |

The per-shell basis gains MORE from data than the default (+0.237 vs +0.182), which is the signature
of real content rather than fitted noise. It has simply not reached the 0.9 threshold at 100k.

Note what the default does at production scale: **99.11% of its power is above 60 A**, and adding
data makes it WORSE (114.7 -> 148.6 A). Its reproducibility advantage is bought by concentrating
further into the low-frequency bulk, so comparing the two on absolute reproducibility flatters it.
Compare deltas, not levels.

**THE DOWNSTREAM METRIC SETTLES IT AGAINST PER-SHELL. Keep it opt-in; do NOT make it the default.**
100k, npreimages=8, spread (mean pair corr of the most-spread triple, LOWER = states differ) and
cmin (worst state's positive mass in its densest 1% of voxels, HIGHER = still a molecule):

| | all 8 delivered | | matched rank 4 | |
|---|---|---|---|---|
| | spread | cmin | spread | cmin |
| default | **0.3583** | 0.0734 | **0.5461** | 0.0734 |
| per-shell | 0.4300 | **0.0854** | 0.7172 | **0.0854** |

Default wins state separation at both ranks by a wide margin. Per-shell wins cmin — its states are
MORE concentrated, so it is not producing blobs — but they are far more correlated with each other,
and its occupancies collapse to ~3-4 substantive states (206/38978/321/3756/256/36482/19817/184)
against default's ~6-7. Match rank before comparing: `spread` is a MINIMUM over C(n,3) triples and
falls with n on combinatorics alone.

**Correction to the scale-N reading above.** Reproducing better with more data shows the content is
SYSTEMATIC, not that it is CONFORMATIONAL. Both halfsets share the same imaging physics, so a
resolution-dependent CTF or alignment residual reproduces exactly the same way. Cross-halfset
reproducibility separates systematic from random; only a downstream metric separates conformational
from systematic — and here it says the extra 24-60 A content does not help the clustering.

Beware also the small-n trap: at npreimages=2 (3 states) per-shell looked slightly BETTER on both
(spread 0.8913 vs 0.9086, cmin 0.0505 vs 0.0454). One triple is not a measurement. The 8-state
comparison reversed it.
