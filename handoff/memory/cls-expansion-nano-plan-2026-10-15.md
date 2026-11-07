---
name: cls-expansion-nano-plan-2026-10-15
description: "user wants cls_expansion on nanoparticle (SINGLE, Pt liquid-cell) class averages with no lp limits, no cross-validation, no FSC shrinkage; plan given 2026-10-15, nothing implemented yet"
metadata:
  node_type: memory
  type: project
  originSessionId: 8a1493a3-3245-4a15-acd2-8cfa98c0abac
  modified: 2026-10-15T15:57:02.399Z
---

2026-10-15 the user asked for a plan to extend cls_expansion (origin/master: simple_flex_cls_expansion,
simple_flex_cls_kernel, simple_cls_expansion_strategy) to nanoparticle class averages. Their stated constraints:
no lp limits (the 20 A FLEX_CLS_BAND_FIRST_LP cap and FRC-derived bands), no cross-validation (the kernel's
coefficient-split CV), no FSC shrinkage (2F/(1+F) extension, fsc2optlp on the averages).

Plan I proposed (not yet approved or built): rank from the Marchenko-Pastur edge of the noise-whitened member
Gram (noise scale from the bulk, Kritchman-Nadler), noise shape from the background outside mskdiam,
neff-targeted kernel bandwidth, unfiltered kernel means, mskdiam mask instead of the inscribed circle, an
in-plane polish before the fit (tseries 2D search has no continuous rotation), exact rotation tangent, rank
and split floor sized to ~20-frame refine2D_nano classes (default neigs 16 / minpop 25 keep every nano class
whole), one budget group per class, a SINGLE `cls_expansion_nano` wrapper (shared memory) with a typed policy
whose default reproduces master bit-for-bit. Python toy (Pt fcc, [110], smpd 0.358) first.

User's answers (2026-10-15): nano classes hold MORE than 20 frames; every particle image is a SINGLE frame
(independent noise between members); the 2D only clusters them, time-restrained; no dataset question for now.
Written up as doc/implementation_notes/cls_expansion_nano_2026_10_15.md (untracked, in the flex2D checkout).

**Why:** the user's next nano work builds on this; the constraints are theirs, not mine to relax.
**How to apply:** do not reintroduce FRC/CV/lp caps in the nano path; Python toy before Fortran. Related:
[[cls-expansion-stage1-ready-2026-10-12]], [[prototype-in-python-before-fortran]], [[flex-pca-must-be-knob-free]].
