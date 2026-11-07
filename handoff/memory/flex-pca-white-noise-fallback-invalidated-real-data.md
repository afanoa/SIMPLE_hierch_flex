---
name: flex-pca-white-noise-fallback-invalidated-real-data
description: "flex_pca silently fell back to a unit white-noise spectrum on all real data until 2026-08-16; every pre-fix real-data conclusion is provisional"
metadata:
  type: project
---

`load_and_validate_sigma` falls back to a UNIT WHITE spectrum when it cannot find
`sigma2_it_*.star`, and the carry-over helper only searches the project directory, its recorded
working directory, and one level up. A cleaned project nested several levels below the refinement
that wrote the star file never finds it. The run prints
`>>> FLEX_PCA SIGMA: using unit white-noise spectrum after background normalization` — this was
read past for the entire campaign.

Measured violation (per-shell whitened residual variance vs the single global scalar; correct
whitening gives ~1 at every shell): Ribosembly synthetic 1.15–1.30, EMPIAR-10076 real
**3.35–4.33**. The synthetic case controls for un-modelled signal, so the ~3x excess on real data
is genuine noise colouring; the debiasing term under-subtracted 3–4x in exactly the shells that
carry conformational signal.

Fixing it (copying the star file where the carry-over can find it) was a large real improvement:
10076 mean best match to deposited classes +0.223 -> +0.679.

**Why:** treat every real-data flex_pca conclusion recorded before 2026-08-16 as provisional,
including the ones in [[flex-pca-downstream-scoring-protocol]] and the RNA/IgG campaign memories.
**How to apply:** check the SIGMA line in the log before trusting a real-data run, and verify the
star file was actually found rather than assuming carry-over worked. Two knock-on effects: the old
"large states" (4961–7016 kA^3) were whitening artifacts and their disappearance is not a
regression, and auto-K prunes much harder under the corrected model so its 0.45 thresholds need
re-validating.
