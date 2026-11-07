---
name: flex-pca-exactness-verification
description: "How to judge whether a flex_pca change is bit-exact — compare MAGNITUDES in flex_pca_embedding.bin, never just whether it differs"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 84327311-352b-40ee-9426-dbaf205684f4
  modified: 2026-08-03T23:12:14.274Z
---

`flex_pca_embedding.bin` is the only deliverable stored at full double precision (`z`,
`contrast`, `latent_second`, residual energies). The `.txt` tables and the MRC volumes are
float32 or formatted, so they hide last-bit changes the cache exposes. Compare it — but
**compare the magnitude of the difference, not merely whether the bytes differ.**

**Why:** on 2026-08-03 a control-flow addition that is a provable no-op in the shared-memory
path (branches guarded by `flex_pca_is_master()`/`is_worker()`, both false) still perturbed
`z` in 90 % of entries. Measured: median relative **7.5e-16**, max absolute 3.1e-13 against
a |z| scale of 3e+02 — one ulp. Adding or removing branches near these hot loops changes
gfortran's FMA contraction and vectorisation, and that alone moves the last bit.

**The mistake to avoid:** an earlier change (restricting the unconditional
`fpl_out%cmplx_plane = CMPLX_ZERO` in `ensure_latent_projection_plane` to the reallocation
branch, worth ~35 % of the projection stage) was REVERTED after seeing exactly this byte
difference, diagnosed as a stale-memory read. It was not — its statistics were identical to
the no-op case (28891/32000 entries, median 7.502e-16). The optimisation was sound and has
since been restored.

**How to apply:** read the doubles at the cache's `z` offset (12-byte header + 8 bytes of
dims + `4*nptcls` bytes of `pinds`) and report max/median relative and absolute deltas
against the |z| scale. ~1e-16 median is codegen noise; 1e-6 and up is a real logic change.
The pipeline IS bit-reproducible run-to-run at fixed `nthr` and fixed binary, so a
difference between two runs of the SAME binary is always real. The band-limit change
(eigenvalues 1.2e-4) remains a genuine failure of an out-of-band-can't-matter argument;
the memset one does not. See [[flex-pca-optimization-2026-08]].

UPDATE 2026-10-24: bit-exactness is no longer required; the user rule is "numerically close is adequate" (see
[[numerically-close-not-bitwise]]). The magnitude judgement above still sizes a difference.
