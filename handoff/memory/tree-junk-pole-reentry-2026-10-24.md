---
name: tree-junk-pole-reentry-2026-10-24
description: "hetero_tree m8c52 cdisc: a split ending junk unit + ONE group makes the rest an unexplored leaf (O r1, M/MF r12); re-entry fix SIMPLE_TREE_REENTER=yes (uncommitted, m8c58t); amplitude rule = 57 % of toy junk but dropped for CAHRA; EM floor moot on cdisc; how to replay a node mid-tree"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T03:15:03.085Z
---

2026-10-24, for the coordinator "2D to 3D heterogeneous state translation" (worktree sharp-proskuriakova-b68280).

- The gap: a split that sets junk unit(s) aside and leaves ONE genuine group gives K* = 1, so the rest becomes a leaf, unexplored.
  - Seen on m8c52 in all 3 cdisc trees, and in none of the 18 older toy trees:
    - O r1: rest 151,683 (junk 0.23); the tree stopped at 3 leaves.
    - M and MF r12: rest 25,799.
  - O r1's kept pole seeded 282 particles that were 96 % genuine. Its 3,686 junk unit (0.998 junk) formed in stages A/B: the small seeded state drew junk in.
  - M r12's pole was junk already at seeding (0.937), with no low-amplitude component.
- The fix: decide() in simple_hetero_tree.f90, behind the env SIMPLE_TREE_REENTER=yes (test switch, default off).
  - The one group re-enters as a child node, with the rest leaf's rows; the junk unit stays a junk leaf.
  - It needs junk holding particles to be set aside, so each re-entry shrinks the group.
  - Uncommitted. Snapshot m8c58t = HEAD + the switches + lever (c) + the fix; gated.
- Junk rules vs truth at the cdisc root: the amplitude rule flags 23,140 at 0.993 junk (57 % recall, 0.11 % of genuine); empty-only flags 0.
  - That is p27's junk-removal edge, but the rule was dropped at m8c30 for CAHRA (9,016 genuine taken; see [[tree-contrast-and-identity-junk-2026-10-23]]).
  - The coordinator agreed it is not to come back.
- The EM floor (EM_MIN 1000, dropped in e22d3d584) is moot on cdisc: the over-split nodes' groups are 18k+ per half.
  - It binds in 13 small m52M_igm splits.
- Runs:
  - re-entry replays stree_port/reent/O_r1 (54199989) and M_r12 (54199990);
  - MF + fix resume fulltoy_m58MFre_cd1 (54200170), which goes to factorial for validate_tree3 against MF and p27 (leaves AND merged units).
- Replay results (one level; the truth has j1 20k and j2 20k):
  - O r1: the rest re-entered as r11 (151,683). Its fallback seeds split it K* = 2 into 82,877 (junk 0.016) and 68,806 (junk 0.48, holding all its j1 and 93 % of its j2).
    That is a junk / genuine separation the empty-only screen cannot make. Closed is not separated (0.085 vs 0.075). Cost 15.6 core-h.
  - M r12: the junk unit is 99.8 % j1, 87 % of all j1. The rest re-entered as r121 (25,799) and split K* = 2 into 16,220 / 9,579: closed 0.142 vs 0.045, j2 0.093 vs 0.195. Cost 7.3 core-h.
  - MF + fix resume (25 core-h, full end phase): r12_s2 became r121, then r1211 -> r12111 + r12112, plus r1212.
    r12112 is a second closed leaf: 3,043 particles, closed 0.672, recall 0.141, map A B C = 1 0 1 at cc 0.953.
    Closed recall in closed leaves goes 0.358 -> 0.499. Leaves 11 -> 13; merge-up 5 units in both; junk unchanged.
    validate_tree3 (factorial session local_133987a6):
    - the merge-up joins r11212 + r12112 into one closed unit, 9,412 rows at 0.745 / 0.499 (MF 0.779 / 0.358);
    - unit NMI 0.398 -> 0.407, against p27's 0.359 / 0.369; purity 0.720 -> 0.725;
    - junk, NMI and the site column are identical to MF; MF's null splits r111 / r1111 / r11111 remain.
  - The coordinator: don't run a switch-off O resume; the candidate default is MF + fix.
- m52M IgM arm (54186423): 85 leaves (m8c17: 21), 370 core-h at 1.42x its price; left to finish for the merged-unit score.
  - Board test: every split passes near its null, so it is not diagnostic.
  - Beyond-EM re-test: 2.4-10x at depth 0-4, about 2.5x deeper.
  - The max-of-20-shuffles null is about a 5 % test per split, with no multiple-testing control.
  - Genuine-split arm-angle eta^2 (stree_port/split_eta2.py): median about 0.24 at depth 5-8, half under m8c17's 0.24 floor.
    The cuts run along mixed angles, like slicing a continuum. 4 of 6 narrow passes have eta^2 < 0.1.
  - 13 splits have a group under 1,000 per half, so the m8c39 floor removal binds on IgM.
  - Its merge-up gives 45 units at unit NMI 0.381, against m8c17's 21 at 0.422.
- IgM slicing variants, on 13 nodes (stree_port/split_null/igm_varC, igm_varB; snapshot m8c60t has the env SIMPLE_SPLIT_NSHUF and a null mean/sd log line):
  - C, 200 shuffles: drops 4 of the 6 narrow passes.
  - C, knob-free Bonferroni on the beyond-EM re-test (alpha 1/21, m = split steps, z > 3.28): drops all 4 narrow splits with eta^2 < 0.1.
    It keeps r111111 (0.24), the deep 0.16-0.19 splits and the controls. It loses r221111 (0.26).
  - Bonferroni on the board test too is over-strict: its statistic saturates at about z 3.
  - B, flex_min_neff=1000: stops nothing; every node falls to the fflex two-group fallback.
- The re-test-only Bonferroni across other trees (stree_port/split_null/bonf_check):
  - it drops 0 truth-supported splits and 8 of 16 nulls;
  - p27cd1 r21 is lost to the build: its re-test statistic is -0.01 on m8c60t.
  - For live use: m = 2 N_root / minpop - 1, fixed at the root. About 50 shuffles suffice, since it needs the mean and sd, not the max.
  - Real data (split_null/m33_c, ole_bc; stree_port/drop_units.py maps each drop's sides onto the merge-up units):
    - m33 (m 460, z_crit 3.71): 24 of 80 splits drop, and 22 of those are moot (both sides one unit). Two change units: r21212 (u19 | 6,719 of u6) and r111122 (u15 | 5,244 of u1, the 133k unit).
    - OLE (m 52, z_crit 3.12): 7 of 30 splits drop, all moot, so the 14 units are kept.
  - Contrast maps vs off on OLE, same node states: K* agrees at 32 of 37 nodes.
    - Root ARI 0.77, against a seed-1 vs seed-2 root ARI of 0.945, so maps moves about 4x what a reseed does.
    - Seed 2 splits r1 and r2 along other axes (ARI about 0).
- m33 re-entry r12212 -> r122121: PASS at z 7.5, with leaves of 10,450 and 14,141 rows and a dipole difference over 24 % of the mask.
- Trap: an sbatch script ending in `[ -n "${X:-}" ] && cmd` exits 1 when X is unset, so SLURM shows FAILED although every step ran.
- IgM motion promotion, CLOSED 2026-10-24: re-ranking the best deep site candidates at the root (stree_port/promote_igm, climb-style starts).
  - Span ratio 0.43-0.51 against the kept sites' 0.44-0.49.
  - The oracle with truth junk excluded (promote_igm_oracle) adds only +0.01-0.02. a3_th's pole becomes kept.
  - The score saturates at one end: the limit is the score, not junk. The deep 0.65-0.70 was conditional on small levels.
- Direction test split_dirmaps2.py with an even-odd noise null:
  - cdisc pairs m52 vs p27: below-root AUC 0.90-1.00, rho with the truth factor 0.45-0.76;
  - "none" splits also recur far above noise (z about 8);
  - on IgM it can't tell slices from real splits (rho 0.26). Commit note (once the user approves): re-entry gives empty-only a junk-separation path without the amplitude rule.

**Why:** the coordinator weighs re-entry as a default; the replays show whether the rest splits on closed / junk type / A.
**How to apply:** replay a node's decision mid-tree with stree_port/stree_reent.sbatch.
- Work on a copy: rsync without fsite/stA0/stageA/seeds/fflex/leaf_*/part scratch.
- Stamp qsys_partition on the copy's p.simple, n.simple and stageB/b.simple: the compenv picks the partition.
- In the node table, put the node at S_SPLIT=12, plus a dummy node "rz" at step 1 with only node_rz/mean/hetero_step.exit.job. That keeps nact > 0, so final() never runs.
- Pass site_layer=no node_proj=no, kill the controller when every real node is 13 or -1, then score with reent_score.py.
- To resume a whole tree, rewrite the absolute paths in its node, leaf and site tables to the copy; otherwise plain() writes into the source.

Related: [[hetero-tree-4p5m-scaling-2026-10-23]], [[never-delete-only-read-or-move]].
