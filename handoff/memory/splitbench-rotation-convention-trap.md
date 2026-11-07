---
name: splitbench-rotation-convention-trap
description: "TRAP: splitbench3.py's default conv 1,1,1,1,1,1 means transpose=1 (einsum 'qj,nij->nqi'); copying its geometry with 'qj,nji->nqi' gives mis-oriented maps that still look reproducible but carry no conformer information"
metadata:
  type: project
---

2026-10-19: tree_belong.py / tree_simmat.py / tree_flat_estep.py (toy55/) copied splitbench3.py's Fourier-insertion
geometry "with identity conventions" but used `np.einsum('qj,nji->nqi', kk0, R)` (R^T k). The bench's default
`--conv 1,1,1,1,1,1` has transpose=1, i.e. `np.einsum('qj,nij->nqi', kk, R)` (R k). With the wrong form every slice
is inserted at the inverse rotation: leaf maps still separate junk and still give clade-like blocks at >40 A, but a
cross-half E-step over the leaves had NMI(conformer; leaf) ~0.005 and the tested-band self-reproducibility of leaf
deviations was ~0.03 (0.4-0.6 once fixed).

**Why:** the six conv flags are a calibrated set (arm_llr.py); "all ones" is not "no transpose".

**How to apply:** when reusing the bench geometry, copy load_half() verbatim (including the convT branch) and run
the truth ceiling first (tree_flat_estep.py --modes truth: E-step with true-group maps; two_50k gives accuracy
0.94, NMI 0.75 at <=20 A). A near-zero ceiling means a geometry bug, not weak data. See [[klem-split-test-statistic-traps-2026-10-19]].
