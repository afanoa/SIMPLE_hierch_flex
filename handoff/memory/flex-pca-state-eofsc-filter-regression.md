---
name: flex-pca-state-eofsc-filter-regression
description: "The per-state eo-FSC optimal filter (4846782f5, 2026-09-01) is why 10180 state maps regressed vs the 2026-08-27 equal-mass maps: it attenuates 7.7k-particle states 2-3x at 20-13 A; SIMPLE_COV_STATE_EOFILT=0 (7d9d38d1f) restores the consensus-FSC low-pass. GMM AUTO on a continuum also gives neff ~N states = consensus copies"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-07T17:00:00.000Z
---

Found 2026-09-07 when the user compared the v3 spliceosome (10180) state maps with
`~/Downloads/flexpca_equalmass_volumes/10180_PC{1,2}_rank{1..4}_state*_n*.mrc`
(2026-08-27 equal-mass slices along one PC, ~7.6k particles per state at box 320).

Two separate causes, measured with `~/mapcmp.py` / `~/refcmp.py` on the box (radial
power vs the consensus `4_abinitio3D/recvol_state01.mrc`, FSC to consensus, and
mean-subtracted best-match FSC):

1. **Filter.** The reference maps' log10 power relative to the consensus is 0.00-0.12
   at every shell to 13.6 A: they carry the project (consensus) FSC low-pass. The v3
   maps fall 0.5-0.9 log10 between 68 A and 13.6 A relative to the consensus because
   `simple_flex_pca_rec3D.f90` applies `fsc2optlp(eo-FSC of the state)` per state since
   4846782f5 (2026-09-01), and a 7.7k-particle state has FSC 0.5 at ~20 A. Commit
   7d9d38d1f adds `SIMPLE_COV_STATE_EOFILT=0` (project-FSC low-pass, eo resolution
   still logged). Default is still the eo filter; the user's judgement is that the
   consensus-filtered maps are the ones that "look right".
2. **State kernel.** GMM AUTO on the 10180 continuum gives states with neff 19k-42k
   of 47k (frac_in_support ~1.0), so every map is near the consensus (FSC 0.82-0.96
   at 13.6 A) and the six differ only in noise. `state_axis=1|2 npreimages=6` takes the
   equal-mass along-path branch (GMM refit skipped), sextiles of ~7.7k particles,
   neff ~13k, FSC to consensus 0.72-0.76 at 13.6 A. States-only rerun needs the
   PRISTINE project: the "requires a project with an eo split" THROW is really a
   state-count check and fires on a project that already carries written states.

Caveat: the mean-subtracted best-match FSC between the reference maps and ANY of our
sets (GMM, axis1, axis2) is <=0.3 at 68 A and ~0 beyond, so the 08-27 run's PCs are
not the v3 embedding's PC1/PC2 (v3 CV cut keeps 2 of 20 axes). Provenance of the
08-27 run is not on fwl-c142339 (no matching dirs); likely the Aug 26-28 ladder on
branch flex-pca-reproducibility. See [[flex-pca-gmm-undoes-equal-mass-placement]],
[[flex-pca-polish-destroys-island]].
