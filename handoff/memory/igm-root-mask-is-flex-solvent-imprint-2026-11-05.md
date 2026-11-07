---
name: igm-root-mask-is-flex-solvent-imprint-2026-11-05
description: "VERDICT 2026-11-05: IgM root focus mask = whole sphere because the cheap children are a NOISE sort imprinting flex's solvent basis (600 A sphere, no envelope); half_variance offset removal NO-GO (core-only mask, arms lost); flex pcg_mskfile envelope cures the symptom but the latent stays blind to 50 kDa arms; nothing ported"
metadata:
  node_type: memory
  type: project
  originSessionId: 720bcd8b-5831-4be2-afd9-171b2c227ae7
  modified: 2026-11-06T04:00:31.325Z
---

2026-11-05 side session cheap-offset. Note: side-session-notes/local_c1bafd06-1da5-43a7-8173-c2f4899054f5/cheap_offset.md.
Batch: igm_cryolo/tree/fixture_igm/cheap_offset/ (scripts, out, env).
- **Not particle background.** The in-mask DC after SIMPLE's norm_noise explains R2 .02-.05 of the root flex latent.
  The fixture latent holds no arm info (R2 <= .001 over all 10 dims). The null twin reproduces every PC (R2 .95-.99).
  PC1 is the noise in flex's contrast estimate (r -.59).
- **The imprint.** The cheap children differ by PC1's own basis volume (polished basis; r .64 fixture, .80 real).
  - It is a ball filling the 600 A sphere, 98 % solvent. The basis is solved on the whole sphere when pcg_mskfile is
    unset: 92x the molecule.
  - half_variance's scale-only fit keeps its offset (+/-.2 of the support mean, cross-half r 1.000), so the 0.2 core
    is the sphere.
  - Same signature on the synthetic arms root (cheap children) and qtoynull / micnull (flex states). Absent on r131,
    qtoy2 and OLE.
- **NO-GO: map-side normalization.** Offset, DC, affine, high-pass or radial removal in half_variance gives core-only
  masks with arms .00. The arms root loses arm 2, qtoynull goes from 20 to 7 sites, and structural gates do not change.
- **Hold: flex envelope.** pcg_mskfile = the node mean's soft support (an existing key) gives offsets +/-.02, cores
  83-91 % in the molecule and masks of 60-71 %. But arm R2 stays <= .024 and the mask is still core-only. Not ported;
  it would need a full gate battery.

**Why:** the whole-sphere mask is a symptom of a root latent blind to 50 kDa arms. No mask fix creates arm states.

**How to apply:**
- Do not re-propose offset, DC or high-pass normalization of the half_variance maps for IgM.
- Large, cross-half-reproducible cheap-child offsets flag a noise-sorted node.
- Arm detection must come from site-first analyses.

Related: [[tree-2d-subclass-variance-mask-fails-2026-11-01]], [[tree-mask-share-core-no-go-2026-11-02]],
[[flex-pca-solvent-mask-dominates-rna-basis]], [[tree-good-pose-null-fails-halves-and-seeding-2026-10-30]],
[[focus-mask-truncates-mobile-tip-falloff-prior-2026-10-21]].
