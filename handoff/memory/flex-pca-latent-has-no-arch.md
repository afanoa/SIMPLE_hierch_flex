---
name: flex-pca-latent-has-no-arch
description: "The IgG flex_pca latent carries no 1-DOF curvature (arch index 0.04-0.06 vs planted 1.0), so tangent/geodesic PCA has nothing to straighten — but it does carry real non-curve structure"
metadata: 
  node_type: memory
  type: project
  originSessionId: 75235947-f21e-4c4b-ab6c-e947c2f45046
  modified: 2026-08-05T16:03:54.835Z
---

Arch (horseshoe) test run 2026-08-05 on three IgG embeddings — `~/igg_flexpca_probe`
(48 comps), `~/flex_em_sweep/iters03` (10), `~/manifoldem/.../em_sweep/c64_p0` (6),
all 100k particles. Instrument: `scripts/flex_arch_test.py`, results in
`reports/flex_arch_test/`.

**Result: no arch, at all three.** Arch index 0.04–0.06 on a scale where 0 is a
no-manifold null and 1 is a planted single-DOF arch. The principal-curve statistic
(held-out VE of a 1-D curve minus a straight line) sits AT or BELOW the null — for
the 48-comp probe run the real value is below null at z = -3.50.

**AMPLITUDE CAVEAT — the original power claim was too broad.** The planted arch
was scaled to the OBSERVED latent variance (~86 % of it), but per
[[flex-pca-latent-is-997-percent-noise]] real conformational signal is only
0.28 % of that variance. So the z = +100..+360 power check demonstrated power
against an arch ~300x stronger than anything real. The correct, narrower reading:
*the observed latent carries no curvature beyond a covariance-matched Gaussian at
the amplitude of the observed variance* — which is almost entirely noise. Whether
an arch exists at the true 0.28 % amplitude was NOT established by that run.

The null arm is unaffected and the negative is not empty: null total covariance
equals the data's exactly, however the signal/noise split is misallocated, so the
"no more curved than a matched Gaussian" comparison stands.

Note the related scale trap: the cached matrix is `G`, and flex_pca's posterior
covariance is `sig2_eff * G^-1`. `flex_arch_test.py` inverted `G` without the
`sig2_eff` factor, which by accident landed ~7-19x under the true within-frame
noise instead of ~1e6 under — which is the only reason the run produced anything.

**Consequence (still standing, with the caveat):** a tangent-space / geodesic /
PGA reparameterization has no large-amplitude curvature to straighten. Same
direction as [[flex-deformation-model-falsified-on-igg]].

**Superseded instrument:** `scripts/flex_latent_structure_test.py` redoes this
with GT-frame ANOVA for the noise scale, runs in both the eigenvalue subspace and
the between-frame (signal) subspace, and plants controls on an amplitude LADDER so
the deliverable is a detection threshold comparable to the measured signal
fraction. GT frame labels need no external file: `frame = (pind - 1) // 1000`.

**But the latent is NOT Gaussian.** The test found real structure of a different
shape: mid/high PCs depend on PC1 through HIGH-degree polynomials while the low
degrees are empty — e.g. iters03 has incremental R^2 of 0.43 at j=7 and 0.64 at
j=8 against a null of 0.03–0.04, with j=2 at only 0.04. Verified not to be a
Vandermonde conditioning artifact: refitting with a Legendre basis on standardized
x gives the same or LARGER R^2 (cond drops from 1e25 to ~25).

A smooth 1-D curve puts its energy in LOW degrees and decays (the planted arm goes
0.96 at j=2 down to 0.09 at j=10); the real data does the opposite. Combined with
the flat principal-curve result — which would also catch a folded or closed curve —
the structure is not one-dimensional. Step/multi-modal or genuinely
multi-dimensional structure is the natural next hypothesis; it has NOT been tested.

Also worth noting from the same run: the probe embedding's participation ratio is
25.2 (48 comps) with PC1 holding only 28% of the top-12 variance, and the signal
spectrum is nearly flat (2728 -> 486 over 12 comps against a noise floor of
188 -> 79). A very diffuse, high-dimensional latent.

Method note worth keeping: at n=100k a z-score measures sample size, not effect
size — the first version of this test fired "arch present" on a z of +10 whose
actual excess was 0.007 against a planted effect of 0.18. Always calibrate against
a planted positive control, not just a null. See [[flex-pca-embedding-cache-format]]
for the per-particle precision that makes the matched null possible.
