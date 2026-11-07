---
name: judge-2d-partitions-by-heterogeneity
description: rank 2D classifications by the heterogeneity they resolve (they guide 3D importance sampling), never by class-average resolution
metadata:
  type: feedback
---

2026-10-03, after I wrote that abinitio2D "outperforms at <=8 A on every data set except the 80S": "your assessment ... is
not considering heterogeneity. Performance should be relative to the amount of heterogeneity presented in the class
averages, as this informs importance sampling, and not really resolution of class averages". Earlier the same day: "we need
to lean into why flex2D is better than abinitio2D more".

**Why:** the 2D partition restricts and weights the 3D search of every particle; a class must hold one state of one view.
Class-average resolution follows class population (at equal population flex2D and abinitio2D averages agree equally with an
independent map), so a view-only partition looks sharper only because its big classes mix states.

**How to apply:** lead comparisons with state separation (NMI/purity, information beyond view), within-view state differences
(diag_views.py pairs: same-view classes differing in state), number of distinct resolved classes, certifiable differences
and multi-state 3D purity; present resolution only as a consequence of population, normalised at equal class population.
Related: [[flex2d-findings-2026-09-21]].
