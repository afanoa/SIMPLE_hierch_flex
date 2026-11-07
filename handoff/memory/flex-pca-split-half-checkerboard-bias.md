---
name: flex-pca-split-half-checkerboard-bias
description: "The checkerboard Fourier-lattice split behind the reliability prior's rho is a half-box shift in disguise; it drove rho to 0 on every paired basis and collapsed the paired-recipe embedding 1000x (fixed 2026-09-07 by a hashed split, SIMPLE_COV_SPLIT_CHECKER=1 restores)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-07T04:05:30.989Z
---

Found 2026-09-06/07 while rerunning the paired half-set recipe on 10028/10180 at
a86bdd2d2 (auto4-restore). Symptom: paired-recipe runs delivered latents with sd 0.07
(auto4 single fit: 11.6), a UMAP ball, GMM AUTO seeing one macro-cluster; every
paired run on the box (shared-memory and distributed) logged `rho= 0.0000` for every
component at the all-N embedding, every single-fit run logged rho 0.24-0.31.

**Mechanism (measured, not a paired bug):** `cov_herm_inner(lhs, rhs, half)` split the
plane lattice by parity of ih+ik (checkerboard). half1 - half2 of conj(U)*y equals
<U, y shifted circularly by (N/2,N/2)>, which is large and pose-dependent because the
particle images are NOT masked (`cov_image_mask_radius` returns 0 by default) and the
basis mask (320 A at box_crop 64) exceeds N/(2*sqrt2). It enters the two half-data
solves with opposite signs -> raw split-half correlations of -0.13..-0.24 on the
paired polished basis, on the data-free init basis, and on a paired basis embedded
inside a single-fit driver. A 6-iteration single fit had enough signal to net +0.1..0.3.
The reliability prior floors rho at 0 -> prior precision = 1/(1e-3 * eig) -> 1000x
over-shrinkage after the re-solve. Not thread-count dependent (bit-identical at 8/16/48).

**Fix:** hashed lattice assignment (balanced 0.5007, autocorr < 0.03) in
`cov_half_parity`, submodule-local flag set in `estimate_covariance_mean`; commit
755500d4e on box branch `auto4-rhofix` (~/SIMPLE-h0901-int2). `SIMPLE_COV_SPLIT_CHECKER=1`
restores the old split. The polar half-set split (even/odd angle index, polar module)
has the same design and is untested; the polar E-step is env-gated off.

**How to apply:** any rho/reliability number logged before 2026-09-07 from the
checkerboard split is biased low; re-derive anything calibrated on it (the RHO_FLOOR,
component ordering by rho, the "rho drives state-target placement" logic). The
auto4 occupancy gate (7881/18064/...) will NOT reproduce bit-for-bit after the fix
because the single-fit rho values change too; a new reference must be measured.
The 20k-window probe harness is `~/run_zhalf_probe*.sh` + `~/zhalf_analyze.py` on the
box (dump via SIMPLE_COV_ZHALF=1 -> flex_pca_zhalf.bin).
