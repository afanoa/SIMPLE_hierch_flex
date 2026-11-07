---
name: flex-pca-state-reproducibility-gap-test
description: "The decisive test for whether flex_pca states are real: same-state vs different-state correlation ACROSS halfsets — the gap, not the raw agreement"
metadata:
  type: feedback
---

To decide whether delivered states carry real conformational signal, do NOT use raw pairwise map
correlation (it measures the shared consensus - the same bug as the merge gate, see
[[flex-pca-merge-fsc-needs-mean-subtraction]]) and do NOT use deviation energy (confounded by
[[flex-pca-states-stratify-by-neff]]).

Use the **cross-halfset gap**: with the per-state even/odd half maps, compare
`corr(even_i, odd_i)` (same state) against `corr(even_i, odd_j), i != j` (different states).
Real states make a state's own half match better than another state's.

Measured on EMPIAR-10049 (RAG SEC), where the answer is a clean null:

| arm | same-state | diff-state | GAP |
|---|---|---|---|
| neigs=4 npreimages=6, GMM on | +0.4240 | +0.4056 | **+0.018** |
| neigs=4 npreimages=6, GMM off | +0.0617 | +0.0611 | **+0.0005** |

The same-state figure alone (0.28-0.59) looks encouraging and is MISLEADING — it is RAG itself,
shared by every state. Only the gap is informative.

**Why:** it converts "absence of evidence" into a controlled negative — the metric demonstrably
sees signal, it just is not state-specific.
**How to apply:** also check the deviation-pair median against the structural floor -1/(K-1)
(mean-subtraction forces it: -0.20 for K=6, -0.067 for K=16). A median near that floor means no
structure beyond the constraint. Difference-of-difference reproducibility across halves is the
same test in another form (measured +0.0099 to +0.022 here, i.e. noise).

**Limit found 2026-10-04 ([[pose-state-coupling-toy-verdict]]):**
- On a homogeneous toy, refine3D_states' four "states" pass this test almost like real ones: gap +1.28 /
  +0.97 / +0.55 at 60-20 / 20-12 / 12-8 A.
- They are view/defocus splits: NMI(state, 50 view cones) 0.163, permutation null 0.0008.
- The gap shows a difference reproduces, not that it is conformational. Pair it with NMI(state, view cones)
  and NMI(state, defocus) against a permutation null.
- **flex_pca states on the same null** (GMM placement at consensus poses): they pass at 60-20 A (gap +1.15 vs
  real +1.29) but not finer (20-12 A +0.26 vs +1.04; 12-8 A +0.04 vs +0.61). Read the flex gap test at 20 A
  and finer only.
- **RNase P refine3D_states final labels** (final_p3d.txt): view NMI 0.027, defocus NMI 0.003, so not view or
  defocus splits.
