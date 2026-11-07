---
name: pole-census-synthetic-validation-2026-10-22
description: "fixed-point census of the 2-class site EM (20 starts) on synthetic truth: discrete states are attractors (rare from its tails, both sides of a switch, a middle cluster); continua give their ENDS (1D: 2, 2D: 3+ edges) -> attractor count is not K; the middle decides; run EMs to convergence (30 it)"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-22T16:37:25.747Z
---

`census.py` / `census.sbatch` / `swing_census.sbatch` (worktree laughing-carson-4e8ca1, routes/toy55); outputs rnasep_klem/factorial/census_*;
results in `doc/implementation_notes/factorial_site_layer_results_2026_10_21.md` ("Fixed-point census on synthetic truth").
- Starts: tails of site-flex z1-z6, middle bands of z1-z3, 3 disjoint random, small tails; group by Spearman; attractor = >= 2 starts.
- Discrete: simsites closed 8 % = attractor from 10 tail starts (98 % pure), random starts never reach it; 70/30 domain: BOTH sides
  attractors (no minority-vs-rest collapse with enough starts); swing three clusters: all three incl. the middle.
- Continua: swing 1D continuum -> exactly 2 attractors (its ends, 6 / 23 deg of 0-30); IgM 2D arms -> 3 edge attractors each.
  Middle starts slide to an end in a continuum but stay in a discrete middle state -> the MIDDLE behaviour separates them.
- TRAPS: "disjoint" must be relative to chance (15 % starts overlap ~15 % at random); Spearman over all rows splits big classes
  before convergence; hold-out filter passes continuum edges (0.72-0.87); EMs NOT converged at 10 iterations (2-3 %/it argmax
  changes) -> 30 iterations: same states, 19/20 starts on them (vs 15).
**How to apply:** count states from attractors only together with the middle-start test, and run site EMs to convergence.
Related: [[factorial-site-layer-test-2026-10-21]], [[landscape-must-be-reproducible-and-splits-data-derived]].
