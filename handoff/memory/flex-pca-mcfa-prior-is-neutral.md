---
name: flex-pca-mcfa-prior-is-neutral
description: "2026-09-18 clean A/B on the 4-dataset suite: the MCFA mixture latent prior (K=16) vs plain PPCA (K=1) is inside run-to-run noise on basis reproducibility, K, occupancies and state resolution on all four molecules"
metadata:
  type: project
---

Same binary (~/SIMPLE_unified_p5, 4c7f4bb7a), same recipe of record, only SIMPLE_COV_EM_MIX differs
(16 = baked default vs 1 = the code's exact plain-PPCA reduction). Jobs 53236059-62 (mix16) and
53236063-66 (mix1); runs under <dataset>/flex_suite/{mix16,mix1}; scorer `suite_score.py unified6 mix16 mix1`.

| dataset | arm | comps >=0.7 / >=0.85 | meanmin | K | states (eo-FSC A) |
|---|---|---|---|---|---|
| ribo   | mix16 / mix1 | 4/4 vs 4/4 | 0.560 vs 0.562 | 4 vs 4 | 8 @ 4.2-4.9 both, occupancies within 3% |
| splice | mix16 / mix1 | 5/2 vs 5/2 | 0.642 vs 0.588 | 5 vs 6 | 5 @ 8.9-9.7 both |
| trpm4  | mix16 / mix1 | 9/6 vs 8/5 | 0.848 vs 0.841 | 9 vs 8 | 5 @ 7.8-8.7 both |
| cnga1  | mix16 / mix1 | 8/3 vs 6/4 | 0.799 vs 0.799 | 3 vs 3 | 5 @ 4.7-6.4, occupancies identical to 3 digits |

The top components (agreement >0.9) are identical between arms on every dataset; the differences are
in the 0.4-0.7 tail, which is the run-to-run band ([[flex-pca-distribution-changes-results]]). The
mixture engaged (log: `MIX distributed reduce: kmix=16`). Wall identical (22-37 min).

**Why:** the earlier evidence for MCFA was bundled (10076: prior + direct-K delivery) or a single
20k row; this is the first isolated test, and it says the prior does not change what is delivered.
The state stage is dominated by the deconvolution + GMM_AUTO placement downstream, not by the
latent prior inside the EM. K=16 is therefore not a differentiator to lean on in a paper without
a labelled-dataset result (10076 with the ceiling instrument) that isolates it.

**How to apply:** keep K=16 as default only if a 10076/Ribosembly labelled A/B shows a gap; otherwise
K=1 is simpler and equivalent. Related: [[flex-pca-mcfa-ab-suite-2026-09-18]],
[[flex-pca-supervised-ceiling-instrument]], [[flex-pca-hierarchical-gmm-auto]].
