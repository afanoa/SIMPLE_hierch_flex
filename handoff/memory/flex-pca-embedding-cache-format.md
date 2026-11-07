---
name: flex-pca-embedding-cache-format
description: "Binary layout of flex_pca_embedding.bin, and the fact that it already stores a per-particle latent PRECISION matrix"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 9671e865-e13f-4298-838f-d873182e5df9
  modified: 2026-08-04T14:17:40.320Z
---

`flex_pca_embedding.bin` is written by every flex_pca run and replayed with `infile=<path>`,
which skips the covariance fit entirely — so it is already the embedding cache; no separate
caching layer is needed for downstream clustering experiments.

Layout from `write_embedding_cache`, `src/main/flex/simple_flex_pca_model.f90:466`.
Fortran `access='stream'`, so **no record markers** — raw bytes, column-major:

```
character(8) 'SIMPLFXC' ; integer version = 3      <- 12-byte header
integer      nptcls, ncomp
integer      pinds(nptcls)                          <- 1-based project particle indices
real(dp)     z(nptcls, ncomp)                       <- the embedding
real(dp)     eigvals(ncomp)
real(dp)     contrast(nptcls)
real(dp)     resid_energy(nptcls), resid_mean_energy(nptcls)
real(dp)     precision(ncomp, ncomp, nptcls)
real(dp)     sig2_eff
```

**The most valuable and least obvious field is `precision`: the per-particle precision
(inverse covariance) matrix of that particle's latent.** It means the estimation noise on
each `z_i` is *known*, not assumed, so heterogeneity analysis can be posed as a
measurement-error model — `z_i | k ~ N(mu_k, Sigma_k + Lambda_i^-1)` — where `Sigma_k` is
the true within-state spread. That separates real conformational spread from estimation
noise without a separate heuristic, and gives the noise floor an eigenvalue test needs.

`precision` is ~99% of the file (`8*ncomp^2*nptcls` bytes), so memory-map it rather than
copying: 1.9 GB at ncomp=50, 7.5 GB at ncomp=100. This, not compute, is what bounds `neigs`.

Reader + `.npz` caching: `scripts/flex_pca_cache.py`. Sanity check of a fresh embedding
against known labels: `scripts/flex_pca_state_check.py`.

`read_embedding_cache` hard-fails on version mismatch, on a particle-count mismatch, and if
any `pinds` entry differs — so a cache is only reusable against the identical particle
selection. See [[flex-pca-exactness-verification]] for how to judge whether two runs' `z`
differ meaningfully (compare magnitudes; ~1e-16 is codegen noise).

**2026-09-16 (SIMPLE_pcg worktree, uncommitted):** the cache is ONE file. After the deconvolution the run
appends a trailing block to `flex_pca_embedding.bin` (magic 'SIMPLFXD', nptcls, ncomp, deconvolved z(n,k) dp,
precision(k,k,n) dp, labels(n) int, noise_scale dp), positioned at `raw_block_end(nptcls,ncomp)`; the raw
reader ignores it. A states-only resume (`infile=<run>/flex_pca_embedding.bin`) adopts the block and skips
the sigma-independent K ladder; older caches fall back to `flex_pca_embedding_deconv.bin` (no longer written).
`scripts/flex_pca_verify/fixedlabel_score.py` reads either.
