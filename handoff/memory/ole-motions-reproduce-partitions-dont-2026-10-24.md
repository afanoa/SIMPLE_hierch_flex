---
name: ole-motions-reproduce-partitions-dont-2026-10-24
description: "OLE 12707 seed 1 vs seed 2 — root-promoted site motions reproduce (|rho| 0.86-0.99, P5-P1 maps 0.92-0.99 at 15 A) while depth-2 partitions are at chance; TRAP root flex identical across seeds"
metadata:
  node_type: memory
  type: project
  originSessionId: 5f12aee4-62b6-4518-9284-7a4df90b8e47
  modified: 2026-10-25T21:12:44.470Z
---

2026-10-24, coordinator follow-up to [[final-estep-ole-real-data-2026-10-24]]. Seeds are 6_hetero_tree (s1) and 7_hetero_tree (s2). Tools: ole_12707/scripts/motion_match.{py,sbatch} and motion_diffmaps.py; output in ole_12707/motion_match/. Cost: 2.1 core-h.
- Sites: each seed-1 factorial kept site was compared with all seed-2 site posteriors (root fsite_* plus node-level site/ down to depth 2) by best |Spearman|.
  - r1_s4 ~ s2 r2_s5: 0.961. r2_s4 ~ s2 r1_s3: 0.988. r_s2 ~ s2 r1_s2: 0.863. The next-best partners are 0.75-0.85.
  - The discovery branch does NOT reproduce: seed 1 found under r1 what seed 2 found under r2, and both promoted it to the root.
  - The node-only site s1 r1_s2 reproduces less (0.73, at s2 node r11).
  - s2 r_s3/r_s4 match an axis that seed 1 rejected at the pole rule (0.95). s2 r_s5, a 5 % pole, has no partner.
- End maps: the seed-2 partner's quintile series reconstructed from seed 1's project with the motion-series recipe ([[motion-series-prototype-2026-10-24]]).
  - P5-P1 correlation at 15 A: 0.97 / 0.99 / 0.92.
  - Null between different motions: about 0.05. r1_s4 and r2_s4 share a component (about -0.7, the same within seed 1).
  - This uses the same images and poses, so it confirms direction; the independent evidence is the rho.
- TRAP: the root flex latent is identical across seeds (rho 1.00 for all components): the seed does not enter it. Cross-seed latent reproducibility must use fflex (per component 0.69-0.99, CCA 0.99-0.80).

- Node level (node_match3/; tools node_match{2,3}.py; 2.3 core-h):
  - Root-inherited node axes reproduce at every size (0.89-1.00).
  - NOVEL node axes: median 0.73 (18-33 % >= 0.8) at nodes of 8k particles or more; 0.67 (0 %) at 5-8k; about 0 below 5k.
  - About half the deep axes have no assessable partner at all.
  - Axes the pole rule kept: median 0.89.
  - End maps for depth 1-3 pairs: 0.90-0.98 at 15 A.
  - So there is a size floor of about 8-10k particles (size and depth are confounded).
- Single-seed gate (2026-10-25):
  - The site layer's existing halves statistics (site_halves self, rel, holdz, polez) do NOT predict cross-seed |r| (AUC 0.51-0.68): they describe the start endpoints, not the converged axis.
  - A real refit does. Rerun the node's main site_rank on random 50/50 particle halves, take each start's high / low groups per half, reconstruct 2 states, score corr(D_A, D_B) at 15 A. On OLE this gives AUC 0.83 for |r| >= 0.8 (Spearman 0.67); 48 core-h (ole_12707/halves_refit/).
  - Fake-direction null (z block shuffled): median |corr| 0.064, MAD-sd 0.076. Its MAX is not a noise level: fake EMs at big nodes converge on real motions (0.38-0.49).
  - Cut = null median + 3 MAD-sd = 0.29: OLE recall 1.00, precision 0.54 (fixed 0.6: 0.86 / 0.60).
  - RNase P m33, 18 nodes of 13-47k at depths 2-7 (rnasep_klem/split_gate_m33/, 110 core-h): the score rises with size at fixed depth (partial +0.48) and falls with depth at fixed size (-0.35). Every node of 13k or more has a passing axis.
  - Gate cost: about 5 core-h per m33 node, dominated by 2-state reconstructions at box 160; m33 groups of about 1.3k particles give cFAR about 0 (signal, not a bug).
  - OLE r2 z1hi: passes the gate, fails cross-seed (0.55) = a site-discovery miss. Seed 2's fflex explains it equally well, but site_dirs=1,3 never seeded it.
- Wider site seeding (site_dirs = all 10 retained fflex components; ole_12707/wide_site/, 194.5 core-h): NO.
  - Site EM converges to 2-3 attractors per mask whatever the start; the 20-start run costs 2.1-4.1x the 4-start one.
  - The r2 miss is site-REGION ranking: the seeds' main masks have Dice 0.
  - In matching regions, wide seeding recovers 1 of 3 missed axes at |r| >= 0.8.
  - Fix chosen instead: every-site ranking (node_motion); infile2 seeding from known axes is the targeted fallback.
  - flex_pca components are NOT eigenvalue-sorted, yet site_dirs=1,3 assumes leading components (logged for the user; flex_pca is the user's code).
- Design trees (2026-10-25: m8c75t / 77t, one motion per split, gated candidates; ole_12707/10_ / 11_tree_design_s*; scorer scripts/tree_pair_verify.py):
  - Cross-seed unit ARI: design 0.319 (10 x 10), same-build baseline 12 / 13 0.254 (10 x 12), old 6 / 7 0.182 (14 x 11). The build and the design each add about half.
    - The design raises matched share (0.38 -> 0.51) more than NMI (0.43 -> 0.40).
    - Leaves alone (baseline 0.138 vs old 0.124) misled me into crediting only the design: compare UNITS.
  - Splits: depth 1 the same in both seeds (0.91 / 0.88); the ungated root differs (-0.71); one private split per seed.
  - s1 r11/c2 is a PROMOTION miss in s2: s2's site posteriors have it at 0.81-0.84, its candidates at 0.29.
  - The gate floor is unstable from few fakes (s1 0.269 vs s2 0.742): fakes converging on real axes inflate the MAD. Fix sent: drop fakes whose |rho| to any start / candidate exceeds 3 / sqrt(N).
- Map-level check vs the user's 31-leaf tree (scripts/map_compare.py on final_maps/raw at 15 A; match = share of the half-map noise ceiling >= 0.7):
  - Only 15 of 6's 31 states reproduce in its own partner 7.
  - The design keeps 10 of those 15; the baseline keeps 13.
  - The design's own states are real (9 / 10 reproduce, 95 % of particles; baseline 9 / 24) and better resolved (median 7.1 vs 10.1 A).
  - The design falls short on the r1_s4 extreme (top leaf mean 0.63-0.66 vs 0.73-0.76) and misses r1212 / r21211 (5.9k / 3.3k).
  - Judge trees by reproducible states + span + resolution, not by partition ARI.
- TRAP: deep-node site posteriors are near-binary. Spearman then gives tie artifacts (exact 1.000 / 0.866) and quintile bins collapse; score a hard axis by rank-biserial (2 AUC - 1) of the other score over its classes.

**Why:** this is the user's real question. OLE partitions below the root do not reproduce (depth-2 at chance), but its motions do.

**How to apply:**
- Describe heterogeneity as per-node motion scores (series and their extremes) rather than leaf slices.
- When testing reproducibility, match motions by score correlation across seeds, not by tree position or site name.
- Never cite root flex agreement as evidence.

Related: [[discrete-splits-capture-the-span-states-that-matter]], [[landscape-must-be-reproducible-and-splits-data-derived]].
