---
name: rnasep-flex-per-state-plan
description: the agreed per-state flex_pca recipe for RNase P (states 1/2/4, solvent mask, npreimages=4) and the refine3D_states masking caveat
metadata:
  type: project
---

Agreed 2026-09-19 for the RNase P (afm_cryoem, 932,496 ptcls) multi-state run
`/mnt/beegfs/elmlund/afan/rnasep_states/`.

**Run flex per state on 1, 2 and 4 only.** States 3 and 5 are **junk sinks**, not
conformations — the user identified them from the maps. They absorb ~22% of assigned
particles and stay at 8.4 Å while 1/2/4 sit at 5.1-5.6 Å. Do NOT read the
population/resolution coupling as neff stratification here; junk simply does not
average coherently.

Recipe per state, on the **native box-160** maps from `calc_final_rec` (never the
box-120 registration maps):
```
selection  infile=<state N 0/1 flags>          # no state= selector exists on flex_pca
automask   vol1=<state N map> smpd=1.464 automsk=yes amsklp=12 nthr=8
flex_pca   vol1=<state N map> pcg_mskfile=<state N mask> npreimages=4 nparts=10 nthr=8
```
`mskdiam` is NOT needed when `pcg_mskfile` is given — the solvent mask replaces the
spherical support the basis is solved on (`rec_backend` already defaults to `pcg`).
Everything else stays default (`neigs=10`, `n_probe_iters=5`); rank was deliberately
left for later.

**refine3D_states used NO envelope mask.** `automsk` defaults to `no`, so the whole run
was `support=sphere posthoc_mask=none phase_randomization=no`. Every per-state
resolution quoted from it is an unmasked spherical FSC and therefore pessimistic; the
flex solvent mask is the first envelope this data sees.

Judge each per-state basis by the `FLEX_PCA XFSC PAIRED record` line (matched pairs and
|cos| lead), not by eigenvalues — the combined run hit |cos| 0.99 with 10/10 matched.
Splitting three ways drops each fit to ~155-190k particles, which raises the detection
floor.

See [[flex-pca-optimal-iters-depends-on-rank]] before raising `neigs`,
[[flex-pca-solvent-mask-dominates-rna-basis]] for why the mask matters,
and [[cnga1-c4-pose-folding-scrambles-protomers]] for the symmetry blind spot.
