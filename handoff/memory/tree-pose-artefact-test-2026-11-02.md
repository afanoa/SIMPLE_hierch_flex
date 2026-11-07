---
name: tree-pose-artefact-test-2026-11-02
description: "pose-artefact splits measured on all saved trees: the per-split flip-fraction difference (dflip, > 30 deg from the node's input pose, within-strata null) separates micelle/PfCRT pose partitions (0.14-0.41) from every legitimate split (<= 0.022; needed splits with truth <= 0.001); 692k root group r4 fires (0.106) and is a confirmed quality blur (half-split check), the closed lineage is clean; proposed flag dflip >= 0.04 & z >= 5, not implemented"
metadata:
  type: project
---

2026-11-02 (this session, note /Users/atifao/SIMPLE/doc/implementation_notes/tree_pose_artefact_test_2026_11.md, scripts
scripts/hetero_tree_pose_study/, data batch2 hetero_tree_20261028/pose_study/; cost < 1 core-h). 180 pairs over 118 nodes (gate arms
base/sub/df on qtoy2, r131 x2, nulls, OLE; 9 micelle fixture trees; real PfCRT tree_umask_s1/umask2_s1) + 105 pairs of the 692k r17d tree.
- Node-level stage-B flip fraction (> 30 deg from n.simple): legitimate <= 0.021 (OLE root), micelle 0.27-0.37, real PfCRT 0.53-0.56,
  692k root 0.094 (fixposes runs 0 by construction). The POSEGATE value 0.165 sits in the decade-wide gap.
- Split-level: the angle-change KS is NOT specific (needed splits reach z 55; real states re-align differently). dflip is: micelle
  delivered artefact leaves 0.19-0.41 (z 49-73), PfCRT roots 0.14/0.17 (z 51/69), PfCRT depth 1 0.02-0.06, OLE root <= 0.022 (z <= 14.5
  by sheer N), qtoy2 / r131 / nulls <= 0.0025. mic3_umask's split is NOT pose-driven (0.0035): a noise split in a bad node = the guard's.
- 692k: root group r4 (100k) holds 17 % turned particles vs 6.5-10 % (dflip 0.106 vs r1, z 75; r4's own split 0.040; view TV 0.17);
  r4's subtree has the tree's most unstable nodes. Closed lineage r1 -> r13 -> r131 -> r1311: dflip <= 0.005, corner split 0.0005.
**How to apply:** report the node flip fraction and each pair's dflip + z at every split; flag dflip >= 0.04 & z >= 5 (merge as a gate);
the fixposes response to an unstable node is the user's call. Look at r4's leaf maps before the next 692k run. Owner of state_split
(math-reworks line) wires it. Related: [[pfcrt-micelle-fixture-tree-pose-scatter-2026-10-29]], [[hetero-tree-algorithmic-audit-2026-11-01]].

UPDATE (r4 follow-up, note section 7): the 692k root group r4 (100,330, 14 leaves, seeded by route A from central site r_s5 whose
kept pole merged three latent tails z1hi+z3lo+z3hi) is a quality-and-alignment state, not junk and not a conformation: stage-B corr
a quartile below its siblings, 17 % turned > 30 deg (7 % near-flips; turning is a low-corr phenomenon in every group), every map
in its lineage with the arm weak in BOTH positions (closed 0.08-0.22, open 0.41-0.72), closed-reference share 0.0069 vs 0.0126 base.
Half-split check RUN (job 54855816, 0.7 core-h, batch2 pose_study/r4111_halves2, scripts/hetero_tree_pose_study/r4111_halves.sbatch):
leaf r4111_s1 at its stage-B poses, closed/open arm: top-corr half 0.263/0.714, bottom half 0.185/0.606, full leaf 0.216/0.642.
Arm density follows quality in BOTH positions, open:closed stays ~3, neither half becomes a state (r2 open leaves 1.2-1.3, corner
1.00): the 'intermediate' r4 maps are a quality blur of an open-arm population, CONFIRMED. Still open: whether r_s5's pole survives
the subspace test's distinctness test in the joint 692k replay. Both the pose flag and a stage-B-corr quality flag would have
marked the root split that made it.
