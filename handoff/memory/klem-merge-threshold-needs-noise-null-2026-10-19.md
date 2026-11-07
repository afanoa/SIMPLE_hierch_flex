---
name: klem-merge-threshold-needs-noise-null-2026-10-19
description: "cross-tree merge statistic (best 3-shell window 40 A..res of the difference FSC) is NOT calibrated at 0.143: independent noise under the same masks gives 99 % 0.39-0.41; NULLCAL=0.99 (noise null) -> null pair 2 states, two-conformer pair 2 (but purity at base rate). Half maps are NOT exchangeable (cross-half EM)"
metadata:
  type: project
---

merge_cross_tree.py options (2026-10-19, toy55, cross_tree_test.sbatch env + XS suffix): SPECMATCH, DOCK, NUIS (34-map
nuisance projection of each difference: scale, shifts, rotations, their 2nd order, radial bands), NULLCAL=q (threshold =
q-quantile of the statistic on independent Gaussian fields under the pair masks/windows). Per-half state counts:
| rule | null_50k (truth 1) | two_50k (truth 2) |
| 0.143 | 5/5 | 8/8 |
| SPECMATCH / DOCK / NUIS | 7-8 / 7 / 6 | 5-7 / 6-8 / 5-6 |
| NULLCAL 0.99 | 2/2 | 2/2 (purity 0.67-0.69 = base rate: NOT the conformers) |
| NUIS + NULLCAL | 3/3 | 4/4 |
The null keeps one split at difference FSC 0.71 (noise max 0.44); in 3 of 4 halves the smaller state is 89-100 % one
family (families are refined separately), not in the 4th. E-O half-map nulls are contaminated: family assignment comes
from the opposite half's maps, so even/odd conformer mix differs by ~7 % within a family (two_50k null tail 0.86).
**Why:** a significance threshold must come from a null, not the 0.143 resolution convention. **How to apply:** use
NULLCAL with the noise null; never build nulls from E/O halves of klem leaves. Related: [[klem-cross-tree-merge-design-2026-10-19]].
