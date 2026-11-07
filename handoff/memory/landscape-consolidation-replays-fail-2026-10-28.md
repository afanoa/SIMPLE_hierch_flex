---
name: landscape-consolidation-replays-fail-2026-10-28
description: "bottom-up equivalence (0-1 merges), reproducibility merge (re-split null too tight) and HDBSCAN-style stability cut (scale-dependent: 14/2/31 units) do not consolidate saved trees consistently; coarse depth cuts agree only by losing extremes"
metadata:
  type: project
---

2026-10-28, adaptive split/merge session; scripts agg_merge.py, repro_merge.py, stab_cut.py, agg_score.py / unit_score.py,
igm_truth.py; outputs /mnt/beegfs/elmlund/afan/hetero_tree_20261028/route9/{agg,repro,stab}; doc section 10.18.
- Bottom-up map-level equivalence (pair mask, delta = smallest kept effect): 0-1 merges per tree. Leaves are by construction at
  least as different as the weakest kept split. Adding localisation (diffuse = mergeable) gives 2-5 merges.
- Reproducibility merge between independent runs: the Bhattacharyya re-split null is 0.994 +- 0.001, so nearly every pair counts
  as distinct (OLE ref 6 merges, strong 0, IgM 1). A significance null cannot consolidate at these N.
- Stability (excess-of-mass) cut: depends on the effect scale. sqrt(rho) helps qtoy2 F1 (0.616 vs 0.574) but equals the leaves
  elsewhere; SNR_d over-merges (OLE 2 units, IgM extremes 2/20); rho is asymmetric between runs.
- Depth-3/4 cuts agree best across IgM runs but deliver 0/20 arm extremes.
**How to apply:** do not propose another significance- or scale-based consolidation without a reason it escapes these failures;
the knob-free candidates are spent. Related: [[null-guards-do-not-stop-692k-deep-splits-2026-10-28]],
[[null-split-diffuse-real-split-localized-2026-10-28]], [[landscape-must-be-reproducible-and-splits-data-derived]].
