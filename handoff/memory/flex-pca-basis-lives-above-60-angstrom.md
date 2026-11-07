---
name: flex-pca-basis-lives-above-60-angstrom
description: "The flex_pca eigenbasis keeps ~93% of its power beyond 60 A and <1% below 30 A even when the full band to 12 A is available — so low-passing costs nothing and the resolution is lost in the column selector, not at lp"
metadata: 
  node_type: memory
  type: project
  originSessionId: 991559ff-d34f-451c-ad5a-ac1ba95b88c7
  modified: 2026-08-10T16:38:28.977Z
---

Measured on IgG 30k (box_crop=64, extent 384 A so shell k <-> 384/k A), power of the bagged
eigenbasis per resolution band, 2026-08-10:

| run | >60 A | 60-40 | 40-30 | 30-24 | 24-20 | 20-16 | 16-12 |
|---|---|---|---|---|---|---|---|
| lp12 (full band to 12 A) | 93.11% | 1.99% | 4.08% | 0.32% | 0.15% | 0.12% | 0.22% |
| lp20 | 91.90% | 2.08% | 4.46% | 0.27% | 0.15% | 0.39% | 0.72% |
| lp28 | 91.25% | 2.52% | 3.99% | 0.67% | 0.48% | 0.39% | 0.69% |
| lp36 | 96.82% | 1.70% | 0.49% | 0.39% | 0.18% | 0.15% | 0.26% |

**Even with nothing removed, the basis is a >60 A object.** So the natural objection to low-passing
("lp=28 throws away detail that could carry heterogeneity") does not apply at lp<=28 on this data —
that band carries <1% of the basis power and was never being used. It DOES apply at lp=36, where the
40-30 A band collapses 4% -> 0.49%; the cross-halfset CV score drops there too (0.650 -> 0.595), so
the CV criterion is tracking real information loss rather than merely rewarding smoothness. That
coincidence is the strongest validation the criterion has.

**THE `|xi| mean= 31.8` / OUTER-HALF WARNING LINES ARE WORKER NOISE — DO NOT TRUST THEM.** A worker
returns from `estimate_snr_volume` with `snr` all zeros (only the master reduces the parts and
finalises), so its greedy selector ranks nothing, ties at 0 everywhere, and keeps the first
admissible lattice point — in the band corner. It then discards that selection (the COLS round reads
the master's `flex_pca_columns.bin`) but logs it first, ten times per run, with the
`Suspect the SNR proxy` warning. Reading those lines as a real defect cost a wrong diagnosis here.
**Always read `flex_pca_selected_frequencies.txt`, which the master writes** — in the same runs it
shows |xi| 1.0-12.1, mean ~8-9 of kmax=32, i.e. healthy low-frequency columns. Fixed in tree by
returning the SNR-stage worker before the selection.

**Where the band restriction really comes from:** the ranking, not a bug. The SNR proxy subtracts ONE
global floor (median of var/dens over shells 0.75*nyq..nyq, ~5.2e-7 on IgG) from a radial profile
that FALLS 1.6e-6 -> 4.9e-7 and flattens beyond sh~22, so only shells below ~22 retain any excess and
the largest excess is at the lowest shells. Columns land at |xi| 1-12 = 32-384 A, and the basis
concentrates above 60 A. See [[flex-pca-per-shell-snr-unconfines-the-basis]] for what changing that
ranking does.

Read with [[flex-pca-lp-is-inert-at-defaults]] (lp inert at shipped defaults) and
[[flex-pca-halfset-fits-lose-fsc-regularization]] (the lp effect is NOT a regularization artifact —
it is bit-for-bit unchanged after that fix).
