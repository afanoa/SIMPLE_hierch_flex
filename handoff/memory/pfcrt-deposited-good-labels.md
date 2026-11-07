---
name: pfcrt-deposited-good-labels
description: "junk-separation instrument for the from-scratch PfCRT picks: deposited 10330 final particles matched to our 592k picks (14,182 of 16,905, corner convention, median 6 px); particle_sieving (5_) DROPS 35 % of them; per-partition AUC/rec@k scorer"
metadata:
  type: project
---

Script `pfcrt_from_scratch/seed_a2d7/pfcrt_good.py stk.txt <deposited star> 300 good.npy <ptcl2D dumps...>`
(copy in the session scratchpad too). Deposited coords are unbinned (x 1918/3838 to our binned mics); our
xpos/ypos are the box CORNER (centre = xpos + 150); mutual nearest within 20 px -> `good.npy` (one 0/1 per ptcl2D
entry, all 592,544 picks share one order across 4_extract / 5_particle_sieving / 6_ / 7_abinitio2D).
Metric per partition: particle-level ROC AUC of the class good-fraction (oracle class ranking), recall of good
when keeping the best classes up to 5/10/20 % of particles, classes >= 3x base rate.

2026-09-30 baselines: 7_abinitio2D (sieved 343k, 200 cls) AUC 0.870, rec@10% 0.59; 6_abinitio2D (all 592k,
300 cls) 0.807 / 0.52; old cold-start flex2D_hetero (400 cls, pruned to 240k) 0.831 / 0.45.
**The sieve keeps only 9,148 of 14,182 good picks (65 %)** -- particle_sieving discards a third of the particles
the depositors refined to 3.2 A. The old flex2D pruning kept 7,434.

**How to apply:** use this as the "junk" ground truth for any 2D-classification or sieving change on PfCRT
(the published-label analogue of 10076). Related: [[empiar-10330-cs-corrected-optics]],
[[particle-sieving-job-model]], [[flex2d-findings-2026-09-21]].

**Full set (seed_a2d6, 592k, 300 cls, 82 min):** seeded flex2D AUC 0.827 vs 6_abinitio2D 0.807, rec@10% 0.57 vs 0.52.
Good particles inside the best classes holding N particles (abs.out):
N = 10k/20k/35k/60k/100k -> sieve+a2d 2319/3965/5440/6906/7933; sieve+seeded 2923/4462/5980/7258/8079;
all-picks a2d 2351/3793/5528/7453/9174; all-picks seeded 2663/4334/6179/8105/9756. Seeded beats its seed at every N;
the sieved route wins small selections, the all-picks route wins from ~35k (the sieve's lost third).
**600 classes (2026-10-01):** seeded 300 -> 600 (seed6_x600, 110 min) AUC 0.831 rec@10% 0.57 vs abinitio2D at 600
(a2d6_600, 106 min) 0.826 / 0.56 -> expansion adds nothing for junk; a2d catches up with more classes. The seeded
gain on PfCRT junk is at matched class count from the seed (0.827 vs 0.807 at 300, 0.888 vs 0.870 sieved at 200).

**Label-free class junk criteria (2026-10-01, class_feats.py in seed_a2d7/ and the scratchpad):** particle-level AUC of
ranking classes by one feature, oracle = good fraction. PfCRT: class-average peak/std (LOW = good; junk = dense dots)
seeded 0.814 (oracle 0.888), a2d 0.809 (0.870), held-out all-picks 300-cls 0.763 / 0.746 (0.827 / 0.807); eo resolution
only 0.693 / 0.513; fixed 4-feature z-sum no better. 10076 (good = published-labelled, 20 % junk): resolution / class corr
best, 0.768 / 0.767 (oracle 0.812); peak/std not in the top. -> the junk signature is dataset-specific (aggregates vs
noise); no single label-free class criterion transfers.

**Class-average 3D (2026-10-01, ab3dc2/{flex,a2d}_{all,top50}: abinitio3D_cavgs pgrp c1 lp 70->8, scored by dock + FSC vs
EMD-20806 with score_vs_emd.sh):** FSC0.143 vs deposited flex_all 51.8 / a2d_all 31.1 / flex_top50 44.4 / a2d_top50 44.4 A
(FSC0.5 78-104 A) -> every model is essentially wrong; top-50 by peak/std (label-free) does not rescue it. The internal
half-map numbers (17.3 flex, 9.7 a2d) are self-consistency, not truth. Class-average ab initio on micelle-dominated,
junk-heavy PfCRT classes is not a usable comparison; the earlier flex 20.7 A (ab3dc_flex, flex2D_hetero 400 cls) was the
best. A meaningful PfCRT 3D test = particle-level abinitio3D on a selected subset.
