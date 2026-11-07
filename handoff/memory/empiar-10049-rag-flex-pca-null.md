---
name: empiar-10049-rag-flex-pca-null
description: "SUPERSEDED 2026-08-20 — the 7-config RAG 10049 null was a POSE + BACKGROUND artifact; native poses killed the megastate and background deflation flattened the nuisance 26:1 -> 1:2, yielding 6 reproducible diverse states"
metadata:
  type: project
---

**Status: the original null is superseded. Read this whole file before citing "10049 is null".**

Original finding (2026-08-17, Blackwell box): seven configurations all null on the cross-halfset
gap ([[flex-pca-state-reproducibility-gap-test]]) — neigs 4/16/20/32, GMM on/off, state_axis=1.
The dissociation was real: the BASIS reproduced (principal-angle cos >= 0.6 for ~32 components)
but no PARTITION did.

**What the 2026-08-19/20 campaign changed:**
- **Poses were the main confound.** The run had used imported poses; a native ab initio
  (abinitio2D -> abinitio3D C2, NU filter) reached **3.58 A** tight-mask (imported 3.69,
  deposited 3.4) and ELIMINATED the megastate: top state 73k -> 10k particles.
- **Background was the second confound.** `DEFLATE_BG` flattens the nuisance spectrum from
  **26:1 to ~1:2**. Separately, pc001 is a single-signed global AMPLITUDE mode and pc002+ are
  genuine dipolar motion modes (clean red/blue domain couples) — the live instance of variance
  ranking picking the nuisance axis.
- **Reproducible states DO exist under that setup.** Over-provision K=24 -> per-state deviation
  e/o-FSC-area floor (kills ghosts) -> amplitude-deflated farthest-point selection gives 12
  reproducible candidates -> 6 representatives with visibly distinct bipolar rearrangements
  (`flex_em/diversity_select.py`). That is cryoDRGN's manual k=100->6 curation, made quantitative.

**What still holds:** 10049 is a CONTINUUM dataset, confirmed twice — merge contracts nothing
(12->12) and deviation-FSC ~= band (21.5 vs 20 A). Het-resolution genuinely pins at 17.8-19.4 A
(see [[flex-pca-derived-lp-is-two-branch]]), 4-5x coarser than the map, so there are no discrete
compositional classes to find and no GT labels to score against. AUTO_NEIGS still wrong here
(caps to its floor of 4). RECOVAR still omits this dataset.

**Why:** "null" was the right read of the measurement and the wrong read of the molecule — the
pipeline was measuring pose error and solvent, not conformation. Always fix poses and background
BEFORE concluding a dataset has no heterogeneity.

**How to apply:** the bundled consensus is the MIRRORED hand — it IS the `_pproc_mirr` volume, and
any flex run against those poses must use it. Next levers: masked e/o reproducibility over the
motion modes for conformational-axis selection, then trajectory delivery along pc002.
