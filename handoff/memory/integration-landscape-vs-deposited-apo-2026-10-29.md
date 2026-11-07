---
name: integration-landscape-vs-deposited-apo-2026-10-29
description: 692k r17d tree vs 55 deposited apo conformers: arm swing covered to ~0.7 A past both ends (conf43 leaf r3111_s2); but by pairwise DOCKED CC our leaves sit inside the deposited landscape (~70 % of its area), no new extremes
metadata:
  type: project
---
User question 2026-10-29: where does our RNase P landscape sit against the published apo set, and do we capture more extremes?
Run on batch2 route9/landscape: fa_maps10.py copy, --dock and --dock --match, both trees' leaves. Figure and tables are in the
worktree at split_merge_replay/landscape/.

**Arm position** (s1 = the deposited swing axis, 92 % of their arm-centroid variance; core-docked, 10 A):
- The deposited set spans -10.1 .. +7.5 A. The integration tree spans -10.9 .. +8.2 A (amplitude-matched variant: -15.8 .. 10.0 vs
  -15.2 .. 9.8).
- Every deposited conformer has an integration leaf within 3 A in the s1/s2 plane.
- **Open end:** r3111_s2 (14.6k) is conf43's mutual best match (arm CC 0.805) and passes the paired test. So does r1221_s1 (-10.9 A).
  The old tree reached conf43 only through junk-like r4 mixtures, which failed the paired test.
- **Closed end:** r1311_s2 (the corner leaf) sits at +8.2 A, past conf18 / 52. It fails the paired test: it is not a close copy of
  any deposited closed conformer. The old tree's r12112 matched conf18 at 0.82.
- **s2 (second axis):** our leaves spread wider, -2.9 .. +4.6 vs -1.7 .. 2.6 A. These are small leaves and pass the paired test.

**Full arm-map detail** (MDS of arm-region distances, only 36 % of the variance in 2D):
- The integration hull is 0.61 vs the deposited 0.79; 5 of 67 leaves lie outside.
- So the tree adds no arm structures beyond the published variety.

**Pairwise docked CC** (the user's preferred measure): every pair rigid-docked at 10 Å, amplitude-matched, box 64.
Script docked_cc.py; run route9/landscape/docked_cc; 1 h on 16 CPUs.
- Deposited nearest-neighbour CC median 0.978. Our leaves' best CC to any conformer median 0.912 (min 0.746); old tree 0.933.
- Projected into the deposited-only MDS, our leaves fall INSIDE the deposited hull (1 of 67 out, by a hair) and cover about
  70 % of its area. They miss its axis-2 extremes, so there are no new extremes.
- In the joint MDS, 50 of 67 leaves sit off the deposited cloud. That is the pipeline offset (uniformly lower CC to every
  conformer; junk-like r4 farthest), not new conformations.
**Subset-size check** (user: we run about 20 % of their data, a pilot; coverage_vs_count.py):
- The coverage gap (best CC to one of our leaves minus the conformer's own neighbour CC) does NOT track deposited particle count
  (Spearman -0.04, p 0.75). The fewest-particle third is covered as well as the largest.
- It tracks atypicality (rho +0.30, p 0.03). The worst-covered are open conformers, including the LARGE conf48 / 49 (140k / 157k).
- Caveat: EMDB counts are curated and enriched (rare extended conformers were deliberately enriched), not natural populations.
  Only more micrographs (about 4x more apo data on the raw share) would settle a subset effect.
Related: [[full692-tree-vs-deposited-benchmark-2026-10-22]], [[integration-closed-two-leaves-2026-10-29]].
