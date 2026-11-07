---
name: view-imbalance-inflates-between-state-metric
description: TRAP (2026-10-14): the reproducible between-state/window map variance (even/odd correlation of state-minus-mean maps) is inflated by any view imbalance between the groups; a view-dominated ordering scored the highest of the night (0.64/0.65/0.38/0.17, one window at 3.1 A FSC, lambda1 0.75); always compare against a view-balanced (within-class) design and a tilt-only partition
metadata:
  type: feedback
---

Found while reconstructing sliding windows along the pooled lp 10 latent on RNase P (section 16.2 of
`doc/implementation_notes/coupled_2d3d_math_formulation_2026_10_13.md`): the direction whose variance is 79 %
class identity (a view partition) gave the largest "reproducible between-window variance" of any family, above the
tracked states in the coarse bands, a single dominant mode (first eigenvalue 0.75) and a window whose half-maps
agreed to 3.1 A. Reconstructions from view-biased subsets differ reproducibly in BOTH halves (anisotropy follows the
view distribution, shared by even and odd), so the metric, the ladder/spectrum readout and per-window FSC all read a
view partition as strong heterogeneity.

**Why:** the even/odd split removes particle noise, not systematic view-dependent reconstruction differences.

**How to apply:** (1) build windows/states view-balanced: inside each 2D class, order the sub-classes along the
coordinate and give every window the same share of each class's members; (2) always run a within-class random
ordering (null) and a tilt-only partition (pure view effect) alongside; (3) report the between-class share of a
coordinate's variance before reconstructing along it; (4) re-read earlier between-state numbers (tracked states,
E-step states, alternation, per-particle bins) with this caveat. Related:
[[rnasep-subclass-window-trajectory-2026-10-14]], [[flex-pca-downstream-scoring-protocol]].
