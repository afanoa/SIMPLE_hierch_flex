---
name: recovar-supplement-unit-conventions
description: RECOVAR supplement formulas assume a unitary transform; porting their absolute constants into SIMPLE literally breaks three separate stages
metadata: 
  node_type: memory
  type: project
  originSessionId: 205ad3f4-4a2f-4250-8466-860b63608957
  modified: 2026-07-31T23:54:29.351Z
---

Measured 2026-07-31 while making flex_covariance follow PNAS 2419140122 S.A-S.E. The
supplement writes its equations in a UNITARY-transform convention that SIMPLE does not
use. Any RECOVAR formula containing a bare `1` or a bare `sigma^2` carries a hidden
convention constant, and porting it literally fails. Three confirmed instances:

1. **eq. (S.9) noise bias.** Supplement says Lambda_i = I after whitening. In SIMPLE the
   whitened per-Fourier-coefficient noise is `sig2_eff ~ 3.88e-6` at box 128 (stable to 4
   digits across runs -- it is a pure box/padding/FFT-scaling constant, not data-dependent).
   Hard-coding the debias coefficient to 1 over-subtracts by ~5 orders of magnitude and
   collapses the reduced covariance. Keep `Dmat = Sbb - 0.5*sig2_eff*SG`; the 0.5 IS exact
   (real half-plane properness), the sig2_eff is not optional.

2. **Algorithm 2 init `R0 = v v^T`, `v = 1e-2*w/sigma^2`.** `w` ("variance of the mean in
   eq. 2") is not definable in SIMPLE's units -- the mean map carries yet another amplitude
   convention. Using the map's per-shell Fourier power with sigma^2=1 puts v^2 below the
   1e-6 density floor, so P pins at the floor, R pins at ~1e6, the per-shell |Sigma|^2 sums
   underflow the DTINY guard, and EVERY column FSC prints 0.0000. Keep the density-relative
   init (commensurate with H by construction).

3. **eq. (S.14) joint (a,z) contrast MAP is NOT IDENTIFIABLE.** Substituting w = a*z, the
   data term depends only on (w,a) while the prior becomes `w'Gamma^-1 w / a^2`, which
   decreases monotonically in a at fixed data fit. Every particle runs to the grid bound --
   measured: contrast mean 2.000, sd 0.000 over 100k particles at A_GRID_HI=2. The
   supplement's own footnote points at Gamma^-1 = 0 (unregularized ML label), where the mean
   term pins a and the problem IS well-posed, but that is a different estimator.

**Why:** these three cost a build+run cycle each (~25 min per run at 100k/box64). The
pattern is worth checking BEFORE porting any further RECOVAR constant.

**How to apply:** port the STRUCTURE of RECOVAR equations, never their absolute scales.
Anything that must be commensurate with H, sigma^2, or the map amplitude has to be measured
in SIMPLE's own convention. See [[flex-covariance-recovar]] and [[flex-covariance-tau-metric-is-uninformative]].
