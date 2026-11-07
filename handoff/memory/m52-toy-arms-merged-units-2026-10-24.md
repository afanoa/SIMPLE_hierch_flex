---
name: m52-toy-arms-merged-units-2026-10-24
description: "m8c52s toy arms vs p27 (2026-10-24): contrast switch needed (switch-off O stops at 3 leaves); M/MF cdisc add null splits that the merge-up collapses (unit NMI 0.40 > p27 0.36); E-step 0.48; IgM M over-splits (45 units vs 21, NMI 0.381 vs 0.422)"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-24T17:57:37.249Z
---

Scored for the coordinating session with [[toy-tree-validation-toolkit-2026-10-24]]. Write-up:
factorial_site_layer_results_2026_10_21.md (m52 sections).
- cdisc s1: O = m8c52s with the switch off; M = SIMPLE_SPLIT_CONTRAST=maps; MF = M + flex_min_neff=1000.
  - O stops after 2 splits: no closed leaf, A-pure share 0.733.
  - M / MF: clean A split, j1 in a named-junk leaf (1.000 / 0.88), closed leaf 0.78 / 0.36 (p27 0.66 / 0.14), but j2 recall only 0.21
    (p27 0.59-0.69; nothing screens it under empty-only).
  - 6 (M) / 3 (MF) null-type splits, all in the A-present branch. Leaf NMI 0.28 vs 0.36.
- Merged units: the merge-up collapses EVERY genuine null pair in M / MF into one A-present unit, so unit NMI is 0.395 / 0.398
  (p27 0.359 / 0.369). It does not collapse p27cd1's r211.
  - Final E-step (one multi-state refine3D pass) then: unit NMI 0.48, j2 R 0.42-0.45, closed 0.83 / 0.74, A-pure share 0.92 (the
    closed state mixes A).
  - Re-entry fix (m58MFre): closed unit 0.745 / 0.499, unit NMI 0.407.
- IgM s1, M: 89 leaves vs the m8c17 baseline's 21. Leaf NMI unchanged (0.421), eta^2 up on most angles, junk recall 0.72 vs 0.92.
  - The merge-up leaves 45 units (the baseline's merge-up merges none). Unit NMI 0.381 vs 0.422, and some units join junk-rich
    leaves across branches.
  - Site abs(rho): a2_ph better; a3_th and a5_th about 0.15-0.21 worse.
**Why:** the delivered output is the merged units, so judge tree changes there. Null splits are harmless on cdisc, but over-splitting
on continuous IgM arms survives the merge.
**How to apply:** score leaves AND units (units_prep.py), and keep IgM in every package. cdisc alone hides the IgM over-split.
Related: [[split-recurrence-across-seeds-2026-10-24]], [[tree-contrast-and-identity-junk-2026-10-23]].
