---
name: flex-pca-lp-is-inert-at-defaults
description: "flex_pca's lp can only coarsen the band, so the shipped box_crop=64/lp=12 pair makes lp a no-op; lp=20 is the first setting that yields a reproducible basis component on IgG"
metadata: 
  node_type: memory
  type: project
  originSessionId: 991559ff-d34f-451c-ad5a-ac1ba95b88c7
  modified: 2026-08-10T15:50:51.499Z
---

`box_crop` and `lp` are NOT two independent knobs. `smpd_crop = box/box_crop * smpd`, so `box_crop`
alone fixes the crop Nyquist, and `covariance_kfromto` applies `lp` only when
`lp > 2*smpd_crop` — lp can make the band COARSER, never finer. The effective band is
`min(2*smpd_crop, lp)`.

On a 128 px / 3.0 A project (IgG), the shipped defaults `box_crop=64, lp=12` give `smpd_crop=6.0`
and Nyquist exactly 12.0 A, so **lp=12 is inert to the last digit**. Verified in the log:
`band kmax=32` (full crop Nyquist), not the 31 an active lp=12 would give.

Measured on IgG 30k, nparts=5, cross-halfset basis cosine (SIMPLE_COV_BAG=1, scored on
flex_pca_bag{A,B}_pc*.mrc), 2026-08-10:

| config | kmax | mean cos | top1 | n>0.9 | basis res_A |
|---|---|---|---|---|---|
| bc64 lp12 (default) | 32 | 0.5192 | 0.868 | 0 | 114.7 |
| bc64 lp16 | 23 | 0.5216 | 0.869 | 0 | 110.1 |
| bc64 lp20 | 18 | 0.6472 | 0.929 | 1 | 106.5 |

lp=20 is the FIRST config producing any component above the 0.90 reproducibility threshold, and it
is not a smoothing artifact — the basis's power-weighted resolution improves (114.7 -> 106.5 A) at
the same time.

Separately: the `columns sit in the OUTER half of the covariance band` warning fires 10x at EVERY
lp. The SNR column selector always lands at ~99% of whatever kmax it is given (31.8/32 at lp12,
22.7/23 at lp16), so lp moves the edge but cannot fix that warning — the selector is the problem.

**Why:** anyone reading `lp=12` as "estimate the covariance to 12 A" is right only by accident on
this box/pixel size and would be wrong on any other.
**How to apply:** a `report_covariance_band` line was added to say which of the two is binding and
warn when lp is inert. Before tuning lp, check `band kmax=` in the log actually changed.
