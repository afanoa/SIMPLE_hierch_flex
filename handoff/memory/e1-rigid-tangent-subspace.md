---
name: e1-rigid-tangent-subspace
description: "2026-09-11: E1 (rigid-motion tangent generators instead of voxel PCA) prototyped on PfCRT and it SHATTERS the covariance detection floor: a 4 deg half-bundle rocking (1.4 A tips) that flex PCA could not see at all (AUC 0.55) is recovered at AUC 0.89 with the amplitude right to 3.9 vs 4 deg; 8 deg gives 0.98. Tools, conventions and the moment estimator"
metadata: 
  node_type: memory
  type: project
  originSessionId: e5141833-c145-408a-ba48-7c3a3ffe0fbd
  modified: 2026-09-12T00:10:59.845Z
---

**Idea:** replace the ~7000-voxel covariance unknown by a few dozen rigid-motion tangent generators. For a region with
soft weight w_r, translation along e is `G = -w_r (e.grad mu)` (per Angstrom) and rotation about e is
`G = -w_r ((e x (x-c_r)).grad mu)` (per radian). Weight applied AFTER the derivative, which matches
"this piece moves, the rest does not" (and the injection's own `w*(rot(V)-V)` form). Dimensionality drops ~100x,
so the detectable amplitude drops ~3x (floor ~ (dof/N)^(1/4) at fixed amplitude^2).

**Tools (batch2):** `~/pcg_tools/e1/gen_build.py {oracle|seg8 K tau}` builds generators at box 300 and Fourier-crops to
box 150 (smpd 2.07); `oracle` = the injection's own two half-spaces, `seg8` = agnostic k-means regions on the consensus
density with softmax weights (no model, no hand segmentation). `reproject` each generator at every particle pose
(EVEN nspace, pgrp, no shifts applied). `gen_fit.py <stack> <out.npz> <mode> [shuffle]` solves per particle in Fourier
space over a 60-5 A band, storing the INFORMATION FORM (b = A^H W y, N = A^H W A) because translation along each
particle's view axis is an exact null direction of A -- never invert per particle. Column 0 is the consensus
projection, so the per-particle contrast is profiled out by Schur complement. `shuffle` = wrong-pose null.
`e1_analyse.py` solves the exact moment identity `sum_i N_i S N_i = sum_i (b b^T - N_i)` (Kronecker, chunked) and
splits the 2-region case into COMMON (rigid body = pose error) and DIFFERENTIAL (conformational) coordinates by
`c = G[m;d]`, profiling the other block out. NOTE a one-half motion splits EVENLY: `m = d = c1/2`, so the expected
differential variance for an injected theta is `(theta/4)^2`, not `(theta/2)^2`.

**Results on PfCRT (16.9k particles, box 150, 60-5 A band):**

| injected motion | flex PCA AUC | E1 matched-filter AUC | E1 specificity (static half) |
|---|---|---|---|
| 4 deg full (1.4 A tips) | 0.55 | **0.887** | 0.547 |
| 8 deg full (2.8 A) | 0.68 | **0.977** | 0.581 |
| 16 deg protein-only | 0.66 | **0.915** | 0.573 |

Permutation null 0.50-0.51 in every case. LDA over all 12 generators, trained even / tested odd: 0.897 at 4 deg.
Fitted end-to-end rotation 3.86 deg vs 4 injected. The tangent model itself is exact at small amplitude
(|cos| with the injected delta volume 0.989 at 4 deg, 0.902 at 8, only 0.460 at 16 -- 16 deg is far from infinitesimal,
which is why PCA's one success was there and E1's advantage is largest where it matters).

**Real-data status (2026-09-11 night, 2-region oracle split, one analysis convention `e1_general.py`):**
the leading conformational coordinate is a RELATIVE TRANSLATION of the two halves. Variances (signal, A^2):
real 0.476 (0.49 A rms) vs wrong-pose null 0.114 (0.24 A); half-set |cos| 0.842 vs 0.010. With the micelle given its own
free rigid motion (2 extra nuisance regions, 24 generators) the excess SHRINKS: real 0.800 vs null 0.335, and the NULL
now reproduces across halves at 0.732 -- so half-set agreement is NOT a sufficient discriminator once a systematic is
present, and part of the apparent protein signal is micelle leakage. Not a per-micrograph/optics systematic
(intra-micrograph ICC -0.02, permutation null the same) and not a mean offset (GLS mean 0.04 A, centring changes nothing).
The decisive null is `make_null2.py`: consensus projection at each particle's own pose plus the RESIDUAL OF ANOTHER
particle (real noise, real micelle, heterogeneity decorrelated from pose).

**Caveats being closed:** the moment estimator is biased by per-particle heteroscedasticity (a global shell sigma makes
noisy particles look like signal), so `RP` (per-particle whitened residual power) is now stored and used to rescale
b and N; and the wrong-pose null (`shuffle`) measures what the estimator manufactures. Real-data differential modes
before those fixes were ~0.5 A rms relative translation of the two halves, halves agreeing at |cos| 0.87, but
per-particle detectability only 0.4 sigma -- do not report until the null is in.

**How to apply:** this is the dataset-agnostic form of the plan's E1/E2 (`doc/implementation_notes/flex_pca_prior_experiment_plan.md`):
regions from the map itself, region count chosen by the independent-half gate, output = per-region motion coefficients.
See [[pfcrt-injection-floor]] for the floor this is measured against, [[pfcrt-poses-are-the-blocker]] for the pose caveat.
