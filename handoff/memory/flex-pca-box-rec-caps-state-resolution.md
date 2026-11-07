---
name: flex-pca-box-rec-caps-state-resolution
description: "flex_pca state maps already resolve to 12-32 A from a >60 A basis, but box_rec defaults to box_crop so they are Nyquist-capped by the COVARIANCE box — raising box_rec is the lever for resolution, not the basis"
metadata: 
  node_type: memory
  type: project
  originSessionId: 991559ff-d34f-451c-ad5a-ac1ba95b88c7
  modified: 2026-08-10T17:49:52.050Z
---

Measured on IgG 100k, npreimages=8, box_crop=64 (smpd_crop 6.0 A, Nyquist 12 A). Even/odd FSC 0.143
of the delivered `flex_pca_state_*` maps:

- default: 12.0 27.4 29.5 29.5 32.0 17.5 32.0 29.5 A
- per-shell: 12.4 24.0 13.2 20.2 12.8 24.0 24.0 13.2 A

**The delivered states reach 12-32 A from a basis that is 93-99% above 60 A.** So a low-resolution
eigenbasis does NOT cap the delivered resolution — it is an embedding device, and the kernel-weighted
reconstruction puts the resolution back. This is exactly RECOVAR's design (Gilles & Singer, PNAS
2025): regularized low-resolution PCA for the embedding, adaptive kernel regression to recover
resolution in the final volumes. Do not "fix" the basis resolution expecting better maps; see
[[flex-pca-per-shell-snr-unconfines-the-basis]], where forcing the basis to 44 A made the STATES
worse.

**RAISING `box_rec` DOES NOT HELP — TESTED AND NEGATIVE (2026-08-10).** I first read the maps sitting
at 12.0/12.4/12.8/13.2 A as "clipped at the box_rec Nyquist, not converged". That was WRONG. Run at
box_rec=96 (smpd_rec 4.0, Nyquist 8 A) and box_rec=128 (smpd_rec 3.0, Nyquist 6 A), same 100k /
npreimages=8:

| state | population | FSC0.143 @ br64 | @ br96 | @ br128 |
|---|---|---|---|---|
| 1 | **774** | 12.0 | 9.6 | 9.8 |
| 6 | **720** | 17.5 | 18.3 | 18.3 |
| 2,3,4,5,7,8 | 10.5k-21.4k | 27.4-32.0 | **unchanged** | **unchanged** |

**The only maps that moved are the two with the smallest populations.** A 774-particle state cannot
support 9.6 A at SNR 0.01 — that FSC is the prior-dominated signature: both halves are regularised
toward the same prior, so they correlate at every frequency regardless of data. Every substantive
state is completely unchanged. The 27-32 A of the real states is a data/method limit, not a box_rec
limit.

**READ FSC AGAINST OCCUPANCY, ALWAYS.** In these runs the two BEST-resolution states were the two
EMPTIEST. An even/odd FSC on a near-empty flex_pca state measures the shared prior, not resolution.

Controlled spread/cmin (Fourier-cropped to a common 64 box, since a finer grid lowers spread and
raises cmin for free): br64 0.3583/0.0734, br96 0.3443/0.0843, br128 0.3429/0.0834 — a small
residual gain, plausibly finer gridding reducing interpolation error, with br128 no better than
br96. But there is NO measured noise floor for `spread`, so 0.014 is not a claim. box_rec=96 is
cheap (1.26 GB machine-wide at nstates=8) and harmless; it is not a lever for resolution.

Per-shell also gives higher-resolution state maps (12-24 A vs 12-32 A) despite worse state
separation — but by the same occupancy argument that is largely its four near-empty states, so read
[[flex-pca-per-shell-snr-unconfines-the-basis]] with this caveat.
