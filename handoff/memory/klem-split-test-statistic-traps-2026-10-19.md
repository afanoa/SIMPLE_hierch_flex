---
name: klem-split-test-statistic-traps-2026-10-19
description: "cross-tree split test traps found on the data-derived tree: group-size bias (compare equal-N reconstructions), nuisance projection CREATES cross-tree correlation, pairing on the tested band double-dips, summed z-score explodes on random subsets; working set = window stat, NUIS=0, MATCHBAND=low, shuffle threshold"
metadata:
  type: project
---

merge_cross_tree.py / dtree.sbatch on null_50k, two_50k, 55_200k pairs, 2026-10-19:
1. Size bias: random groups of different N (2.9k vs 17k) reproduce a difference at FSC 0.91 across independent trees
   (shuffled-label control). Reconstruct every compared group from the same N (neq_subsample.py).
2. NUIS (34-38 nuisance-derivative maps fitted per tree): the fit coefficients pick up the shared coarse component and
   copy it into every frequency the basis covers -> projected shuffled differences correlate at 0.2-0.3 over shells
   12-27 (raw ~0). Never project; NUIS=0.
3. Pairing the trees' groups on 15 A deviation correlation overlaps the tested band; MATCHBAND=low pairs on shells
   below BANDLO (Fourier, one shell gap).
4. STAT=z (coefficient-weighted FSC sum over the band) gives -30..+30 on random subsets even without NUIS (tiny
   systematic high-res FSC x thousands of coefficients). Use STAT=window (best 3-shell mean, 40 A..res).
5. Threshold = max pair statistic over NSHUF label shuffles at equal N (~0.2-0.45; Gaussian-field null too low).
With 1-5 the null root -> 1 state; two_50k root -> 5 (intermediate groups along the real axis: archetype residual 0.39 vs
threshold 0.37); 55_200k root -> 4 broad children (one 0.79-0.87 in the big reference clade). Baseline: the old chain on
55_200k had leaf purity 0.10, strict recovery 8/55; the per-image ceiling at refined poses saturates at ~70 % even with
1000 pooled images, so at 200k only clade-level states are supported. Related: [[klem-merge-threshold-needs-noise-null-2026-10-19]],
[[klem-forced-family-split-and-offset-mode-2026-10-19]], [[arm-llr-ceiling-needs-fitted-map-scale]].
