---
name: flex-pca-transport-rank-is-not-the-barrier
description: "Ribosome ratcheting needs only rank 2-3 in an optimal linear basis, so the Kolmogorov n-width is NOT what blocks flex_pca on transport-dominated data"
metadata:
  type: reference
---

Measured directly on the EMPIAR-10028 consensus (box 64, smpd 7.54 A, lp 20 A): rotate the
volume over a range of angles, band-limit, SVD the deviations, read the rank needed.

| theta_max | transported mass | rank@90% | rank@95% | rank@99% |
|---|---|---|---|---|
| 1 deg | 6.5 % | 2 | 2 | 2 |
| 5 deg | 27.3 % | 2 | 2 | 2 |
| 10 deg | 45.5 % | 2 | 2 | 3 |
| compositional control (lobe occupancy 0->1) | 6.7 % | 1 | 1 | 1 |

**Ratcheting over +/-10 deg needs rank 2-3, not 30.** A rigid rotation is a *one-parameter*
family, and a smooth 1-D manifold is intrinsically low-rank no matter how much mass moves.
The Kolmogorov n-width barrier described in [[transport-dominated-mor-is-the-right-framing]]
bites for multi-parameter families or uniform approximation over large deformations — it is
**not** what blocks a single ratcheting coordinate.

Consequences:
- "large transport needs many components" is WRONG as an explanation of why flex_pca fails
  on 10028 while working on 10076. Do not use a rank-budget argument for this.
- `neigs=40` is not protective. Forty components cannot rescue a signal that sits below the
  per-particle noise floor, and more components mean more parameters estimated per image.
- The real blocker is variance ORDERING: PCA ranks by variance, and solvent + envelope
  nuisance outranks the ratcheting deviation. See
  [[flex-pca-nuisance-deflation-flags-unset]].

Open hypothesis (NOT yet measured): detectability may scale with **signedness** rather than
transported mass — an occupancy change is single-signed and projects strongly onto
consensus-shaped directions, while transport is dipolar (mass moves A->B, net ~0) and so
carries less variance for the same displaced mass. That would explain the observed
10076 (works) > 10049 (partial) > 10028 (fails) gradient.

Scripts: `~/nwidth.py`, `~/span.py`, `~/latent_rot.py`, `~/score.py` on fwl-c142339.
