---
name: flex-pca-rename-and-cleanup
description: "flex_covariance became flex_pca on 2026-08-03 — CLI, module names, all 27 output filenames, and the lot=/diffeo= parameters changed"
metadata: 
  node_type: memory
  type: project
  originSessionId: 34cf3bd2-fff1-4dd4-bda2-2ec6d52c01a4
  modified: 2026-08-03T22:45:22.360Z
---

Committed 2026-08-03 as `working flex_pca implementation` on master.

**What changed at the interface:**
- Program `flex_covariance` -> `flex_pca`; modules `simple_flex_covariance_{model,columns}` ->
  `simple_flex_pca_{model,columns}`.
- All 27 output literals `flex_covar_*` -> `flex_pca_*`, including the resume cache
  `flex_pca_embedding.bin`. An `infile=` pointing at an old `flex_covar_embedding.bin` will not
  resolve.
- `lot=` and `diffeo=` are **removed** — undeclared and unregistered, so an old command line
  carrying them now fails rather than ignoring them.
- `simple_flex_diffeo.f90` deleted. `simple_flex_lot.f90` keeps only the transport core plus
  `lot_pullback_metric`; the descriptive `lot_analyse_*` half is gone.
- All references to the external reference implementation were stripped from the Fortran
  (comments rephrased to "the reference") pending a proper citation.

**Silent breakage the rename left behind (verified 2026-08-03):** the Python scorers were NOT
updated with the Fortran. `scripts/score_latents.py:9`, `scripts/flex_covariance_validate.py:95`
and `scripts/flex_gt_volpc1.py:157` still name `flex_covar_coordinates.txt`, while flex_pca now
writes `flex_pca_coordinates.txt` (`simple_flex_pca_model.f90:2201`). Two of the three build that
path themselves, so they fail against any current run. Check the scorers before trusting any
"run it and compare" step.

**Defaults are not the validated configuration:** shipped `neigs=6`, `min_neff=2000` and the
clamped `npreimages` differ from every recorded parity result, which used `neigs=20`,
`min_neff=100-300`, `npreimages=20`. Comparisons run at defaults reproduce nothing on record.

**Local-only, deliberately never committed:** the ~14 audit docs under
`doc/for_developers/ideas/` and everything Python under `scripts/`. That means the inline
comments in the flex modules are the ONLY committed record of the measured results, which is why
the comment trim kept every sentence containing a number even while cutting density from ~28% to
~10%.

**Unfinished lead:** `lot_pullback_metric` produces a metric that `build_covariance_state_weights`
already accepts as its optional `zmetric` argument, but nothing wires the two together. Measured
offline (run 49, simplified kernel): mean weighted GT spread per state 4.049 -> 3.367 at blend 0.5.
See [[flex-analysis-diffmap-never-worked]] and [[flex-covariance-basis-defect]].
