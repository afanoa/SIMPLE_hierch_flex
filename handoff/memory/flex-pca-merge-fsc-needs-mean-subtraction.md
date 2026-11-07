---
name: flex-pca-merge-fsc-needs-mean-subtraction
description: "The merge gate correlated RAW state maps, so it measured the shared consensus; subtracting the ensemble mean first makes a single 0.98 threshold transfer across datasets"
metadata:
  type: project
---

Root cause of "the merge threshold is arbitrary": the consolidation gate computed FSC between
**raw** state maps, which are dominated by the density all states share. The statistic measured
consensus, not difference — it pushed every pair to the top of its range and crowded correct
fusions (0.9983–0.9995) against harmful ones (0.9848–0.9928). No fixed threshold could transfer.

**Fix (2026-08-17 handoff): subtract the ensemble mean before the FSC**, with reliabilities from
the same deviations. Ribosembly should-fuse 0.983–0.999 vs must-stay-apart 0.413–0.970;
EMPIAR-10076 0.692–1.030 vs -0.976–0.878. A single gate at **0.98** fuses only true duplicates on
both, zero false fusions. Default moved 0.95 -> 0.98.

Refuted with numbers, do not rebuild: difference reproducibility (two subdivisions of one basin
differ systematically, so the difference reproduces too — 0.730–0.936 vs 0.937–0.999, no margin);
between/within in noise units (AUC 1.000 on Ribosembly, fully overlapping on real data, confounded
by per-state SNR); adaptive largest-gap gate (cuts at 0.787 on Ribosembly, ARI 0.948 -> 0.786).

**Why:** the threshold was arbitrary because the statistic was wrong, not because thresholds are
unknowable.
**How to apply:** with auto-K active the merge should default OFF anyway — auto-K delivers a stable
count and the leftover duplicate pair costs ~0.002 ARI. Keep the pair-ratio distribution as the
regime detector: max and spread separate discrete (Ribosembly 0.999/0.71, 10076 0.993/0.47) from
continuum (IgG 0.963/0.15). See [[flex-pca-white-noise-fallback-invalidated-real-data]].
