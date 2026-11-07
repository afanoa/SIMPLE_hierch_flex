---
name: flex-pca-distribution-changes-results
description: "flex_pca is deterministic in one process (20k harness reproduces d' 2.48 / purity 0.50 / rank 6 across dates and binaries) but distributing it over 10 parts gives 2.21 / 0.48 / rank 7: the part-reduce summation order is amplified by the EM, so equivalence tests must compare arms with the SAME distribution, and the per-component rank cut is the fragile piece (principal angles between the half subspaces are identical in both arms)"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-08T18:44:39.239Z
---

Measured 2026-09-08 on the 20k harness (10028, v3 paired recipe): shared memory 80 threads gives
d' 2.48, kNN purity 0.50, rank cut 20 -> 6, identical to the 09-07 harness on a different binary;
nparts=10 nthr=8 on the same binary gives 2.21 / 0.48 / rank 7. The 400-particle gate shows the same
thing: s4 vs s1 z differ 27%, p4 vs p1 134%, while each arm is bit-identical to its own baseline.

**Why.** Only the summation order of the part reduce differs; the half-set EM amplifies that into
a different basis, and the CV rank cut (matched per-component |cos| >= 0.5 between the two half
bases) flips on borderline components: pairs at 0.94/0.72/0.72 in the distributed arm vs
0.60/0.48 (component 2 matched to B's 20th) in shared memory. The principal-angle cosines between
the two half SUBSPACES, printed one line above in the merge, are nearly identical in both arms
(0.988, 0.982, 0.977 ... 0.894 over 13 dims): the subspace is reproducible, its individual
components are not (near-flat spectrum inside it, components rotate between halves).

**How to apply.** (1) Never judge cache/residency/kernel equivalence across arms with different
nparts; compare to the raw arm with the same distribution (10 x 8 raw at 20k = 2.21/0.48/7).
(2) The rank should be counted from the principal angles with a threshold calibrated by a null
(both halves start from the same iteration-1 mixture, so some agreement is shared initialisation),
and the merged basis taken as the principal vectors -- a contained change in em_pairmerge; the user
agreed the current metric is wrong (2026-09-08). Related: [[flex-pca-speed-levers-2026-09-08]],
[[flex-pca-basis-rank-over-specified]].

**CORRECTION 2026-09-20:** one process is NOT bit-deterministic run to run either, at least on batch2 (gcc 15.2,
-O3 -march=native, OpenMP 16 threads, box 96 toy, 40.5k ptcls): the same binary twice gave pre-EM sig2 2.4650 vs
2.4649 e-1 and final eigenvalues 1141/714/391/377 vs 1129/705/411/392 (1-4%), the same spread as master-vs-branch or
build-vs-build. Threaded reductions. Judge code changes by first-iteration numbers agreeing to 4 digits AND a repeat
run of the reference, never by endpoint eigenvalues alone.

**2026-09-30:** the random sigma2 bootstrap (calc_pspec, 25k cap) was also active on that 40.5k toy; with SIMPLE_SEED an nthr=1 run is bit-identical, so the "threaded reductions" cause above is unproven. See [[flex-pca-ab-needs-simple-seed]].
