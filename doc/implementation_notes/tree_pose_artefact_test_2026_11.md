# hetero_tree: pose-artefact splits, a measurement and a test (2026-11-02)

Question (audit 2026-11-01, gap 9): every node re-aligns its particles locally, so a split can be carried by a systematic pose
difference between the children rather than by structure. The micelle fixture showed about 30 % of a root's particles turning
more than 30 degrees at the first refinement and splits that looked aligned-versus-misaligned; a node-level pose gate prototype
(`SIMPLE_TREE_POSEGATE`, off) made the nulls pass. This note measures the per-split signature on every saved fixture tree, the
micelle fixture, the two real PfCRT trees and the 692k RNase P tree, against a null, and proposes the test.

Read only. Scripts: `scripts/hetero_tree_pose_study/pose_split.py` (the measurement, modelled on the micelle session's
`flip_survey.py`), `pose_split.sbatch`, `analyze_pose.py` (the join to the effect-size table and the tables below). Data on
batch2: `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/pose_study/` (`pose_gates_pairs.tsv`, 180 pairs over 118 nodes;
`pose_gates_nodes.tsv`; `pose_rnasep_*.tsv` for the 692k tree). Cost: one 16-core job of 47 s for the gates and fixtures.

## 1. What is measured

For every node that kept its input project (`n.simple`, the poses it started from) and its stage-B project (`stageB/b.simple`,
the poses and states it refined to), per particle:

- the angle between the input and the stage-B rotation, SIMPLE's `euler2m` convention, minimised over the point group for c2
  trees (an equivalent pose of a symmetric volume is `R S_k`, i.e. the first Euler angle plus 360 k / N);
- the change of the in-plane shift.

For every candidate pair of children (the delivered split groups of `split/split_children.txt` where the node split, the
stage-B states where the split test refused):

- `flip_a`, `flip_b`: the fraction of each child's particles that turned more than 30 degrees; `dflip` = their difference;
- the KS distance between the two children's angle-change distributions and between their shift-change distributions;
- a z for each against 200 shuffles of the child labels inside the split test's own nuisance strata (40 Fibonacci view bins of
  the input pose x defocus tertile), so a view-driven split cannot register as a pose-driven one;
- the total-variation distance between the children's view distributions, as the control.

Labels come from the rare-state note's effect-size table (`detlim2/splits/saved_trees3.tsv`, joined by tree, node and child
sizes): needed / unneeded by the truth (Cramer's V >= 0.3 and permutation z >= 5), the nulls unneeded, r131's root and r1
needed, OLE without truth; the micelle fixture trees are the positive control; the real PfCRT trees are the suspects.

Trees: the gate arms base, sub and df on qtoy2, r131 (both seeds), qtoynull, olenull, micnull, OLE (21 trees); the micelle
fixture's trees mic3 / micnull x nomask / umask / order, plus the three fixposes runs fp_order, gate_order, mixB_order
(9 trees); real PfCRT `tree_umask_s1`, `tree_umask2_s1` (23 and 26 nodes); the 692k recipe17d tree (84 nodes, section 4).

## 2. Node level: the fraction of particles the node's own refinement turns more than 30 degrees

| trees | stage-B flip fraction, max over nodes |
|---|---|
| qtoy2, r131 both seeds, qtoynull, olenull (all arms) | 0.0009 to 0.0020 |
| OLE (all arms, c2-reduced) | 0.021 at the root, 0.0035 median |
| micelle fixtures mic3 / micnull, nomask / umask / order | 0.27 to 0.37 |
| micelle fixposes runs (fp_order, gate_order, mixB_order) | 0.0000 by construction |
| real PfCRT tree_umask_s1 / tree_umask2_s1 | 0.56 / 0.53 at the root, 0.05 / 0.06 median over the tree |

This reproduces the 2026-10-29 flip survey (321 nodes: legitimate stA0 <= 0.025) at stage B and adds the real PfCRT trees. The
gap is more than a decade: no legitimate node reaches 0.03, no micelle or PfCRT root is under 0.27. The prototype's gate value
of 0.165 sits in the middle of it.

## 3. Split level: do the children differ in how far their particles moved?

The angle-change KS statistic is not specific: real conformational splits re-align their two states differently, so needed
splits reach z 55 on it (median 12) and it cannot tell a pose artefact from a state. The specific quantity is the flip
fraction difference `dflip`, the share of particles that left their input pose by more than 30 degrees, compared between the
children.

| set | pairs | dflip | z (within-strata null) |
|---|---|---|---|
| needed splits with truth (qtoy2, r131 root and r1) | 20 | <= 0.0011 | <= 4.8 |
| unneeded splits with truth and the refused null candidates (qtoynull, olenull) | 8 | <= 0.0025 | <= 3.0 |
| OLE, all arms (no truth; its root split is real by the maps) | 34 | <= 0.022 (root pairs 0.011 to 0.022) | up to 14.5 |
| micnull in the gates, refused root candidate (guard-merged) | 3 | 0.071 (flips 0.31 / 0.24) | 22 |
| micelle fixture, delivered splits: mic3_order leaf 3 against leaves 1 and 2 | 2 | 0.41, 0.38 (flips 0.65 against 0.24 / 0.27) | 73, 60 |
| micelle fixture, delivered splits: micnull_umask leaves 2 and 3 against leaf 1 | 2 | 0.19, 0.28 | 49, 52 |
| micelle fixture, delivered: micnull_order, mic3_order leaves 1 against 2 | 2 | 0.049, 0.030 | 15, 10 |
| micelle fixture, delivered but not pose-driven: mic3_umask | 1 | 0.0035 | 0.15 |
| micelle fixture, refused candidates (nomask) | 2 | 0.16, 0.12 | 45, 35 |
| real PfCRT, root split (both trees) | 2 | 0.14, 0.17 (flips 0.49 against 0.64, 0.46 against 0.63) | 51, 69 |
| real PfCRT, depth 1 | 4 | 0.021 to 0.063 | 6 to 18 |
| real PfCRT, depth 2 and below | 43 | <= 0.037, median 0.010 | mostly under 5 |

Readings (the user's caveat first: real PfCRT has no heterogeneity, so every split in its trees is false by assumption; it is a
real-data pose null showing that the statistic fires on a pose-unstable dataset, not a test of discrimination against real states,
which rests on qtoy2, r131 and OLE):

- **The artefact mechanism is real and specific.** In the micelle fixture the tree delivers a leaf that is mostly misaligned
  particles (65 % turned more than 30 degrees against 24 % in its siblings), and the real PfCRT roots split 49 against 64 %
  flipped. No legitimate split comes within a factor of 6 of the weakest delivered artefact (0.19) and the needed splits with
  truth sit two orders of magnitude lower.
- **Not every split in a pose-unstable node is pose-driven.** mic3_umask's split has equal flip fractions (0.32 / 0.32,
  z 0.15); it is a noise split of a bad node, which the locality guard, not a pose test, has to catch. The node-level
  fraction (section 2) says the node is unsafe; the split-level statistic says which of its splits are pose partitions.
- **Depth attenuates it.** Below the PfCRT roots the flip fractions fall to 5 to 16 % and the differences to 0.02 to 0.06 at
  depth 1 and under 0.02 further down: the first split absorbs most of the misalignment, and the deeper splits are no longer
  pose partitions even though the tree is still re-aligning 5 to 10 % of each node.
- **OLE is the only legitimate tree with a measurable signal**: its root children differ by up to 2.2 points in flip
  fraction (1.3 against 3.5 %) at z up to 14. The absolute difference is 6 to 8 times below the weakest artefact, and its
  root split is real by the maps (P / N 30 to 61 in the rare-state note). A z alone would flag it; the absolute difference
  does not.

## 4. The 692k RNase P tree (recipe17d, 84 nodes, 105 pairs; one 16-core job of 2 min 20 s)

| node | particles | stage-B flip fraction | split | dflip | z |
|---|---|---|---|---|---|
| r (root, 4 groups) | 692,392 | 0.094 | r1 against r4 | 0.106 (flips 0.065 against 0.172) | 75 |
| | | | r2 against r4 | 0.085 | 54 |
| | | | r3 against r4 | 0.072 | 34 |
| | | | r1 against r3 | 0.034 | 33 |
| | | | r1 against r2 | 0.022 | 19 |
| r4 | 100,330 | 0.037 | r41 against r42 | 0.040 (0.022 against 0.062) | 30 |
| r31 | 79,973 | 0.014 | | 0.020 | 19 |
| r1, the closed lineage | 252,458 | 0.008 | three groups | <= 0.005 | <= 8.5 |
| r13 | 80,775 | 0.004 | | 0.004 | 8.4 |
| r131 | 44,334 | 0.0015 | | 0.0007 | 1.4 |
| r1311 (the closed corner's node) | 22,464 | 0.0013 | r1311_s1 against r1311_s2 | 0.0005 | 0.2 |
| all other nodes | | <= 0.024 | 97 pairs | median 0.002, 90th percentile 0.009 | |

Readings:

- **The root is moderately pose-unstable** (0.094, four times the legitimate maximum, a third of the micelle fixture), and
  the instability is concentrated: the fourth root group r4 holds 17 % turned particles against 6.5 to 10 % in r1, r2 and r3.
  Against r1 that is a flip difference of 0.106, inside the real PfCRT root range (0.14 to 0.17) and six times above any
  legitimate fixture split. r4's own split repeats it (0.040), and r4's subtree holds the most pose-unstable nodes of the
  whole tree (r42 0.024, r421 0.018, r412 0.016). Its view distribution also differs most from its siblings (TV 0.17).
  Whether r4 is a pose artefact, a junk-rich pool whose particles cannot be aligned, or a real state that re-aligns far,
  the pose statistic cannot say; the rare-state note's effect-size table does not cover the 692k tree. It is the one place
  in the production tree where the proposed flag fires, and it deserves a look at r4's leaf maps before the next production
  run.
- **The closed lineage is clean at every level.** From r1 down to the corner node r1311 the flip differences are 0.005 and
  below, the corner split itself 0.0005 at z 0.2. The closed corner is not a pose partition.
- **The remaining splits** sit where qtoy2, r131 and OLE sit (median 0.002, 90th percentile 0.009); r1 against r3 (0.034)
  and r31 (0.020) are the two nearest the line, both below the proposed 0.04.

## 5. The test

Two statistics, both free (the poses are already in the node's projects; the controller reads `n.simple` and `stageB/b.simple`
at decide time, or `state_split` computes them from the same projects):

1. **Node pose stability**: the fraction of the node's refined particles that turned more than 30 degrees from their input
   pose, symmetry-reduced. Report it at every node. A node above 0.165 (the existing prototype's value, in the middle of the
   decade-wide gap between 0.021 and 0.27) is pose-unstable; the existing `SIMPLE_TREE_POSEGATE` response, switching the
   tree to fixposes, made the nulls pass on the fixture and is the candidate action.
2. **Split pose partition**: for every candidate pair, `dflip` with its z against the within-strata label shuffles. Flag a pair
   when `dflip >= 0.04` and `z >= 5`; as a gate, merge it. On the saved trees that flags every delivered micelle artefact
   split (0.19 to 0.41), the refused micelle candidates (0.07 to 0.16), both PfCRT root splits (0.14, 0.17) and PfCRT
   tree_umask_s1's r2 (0.063), and no split of qtoy2, r131, OLE or the nulls (maximum 0.022, a margin of 1.8x to OLE's
   root and 40x to the needed splits with truth). It leaves the mic3_umask split alone, correctly: that one is the guard's.

What the test does not do: it does not catch a noise split in a pose-unstable node whose children are equally misaligned,
and it does not judge whether the state it protects is worth a state. Those stay with the locality guard and the effect-size
floor. It also cannot say, on its own, whether a flagged split should be merged or the node re-aligned with fixed poses;
on the fixture the second was the response that made the nulls pass, and on real PfCRT it is the one the memory already
recommends (treat the leaves as pose artefacts).

## 6. Verdict

Go for the measurement, as two reported numbers at every split (the node's flip fraction and each pair's `dflip` with its
z), written to the node record beside the split statistics; go for the split-level flag at `dflip >= 0.04, z >= 5` as a
gate once the owner of `state_split` wires it (the math-reworks line owns that file; this note proposes, it does not
implement); the node-level fixposes response stays the user's call, since it changes the tree's refinement rather than a
decision. Cost of the study: under one core-hour.

## 7. Follow-up: what r4 is (2026-11-02, user's question)

Sources: the tree's node and leaf tables, the root's split log and seeding log, the saved closed-node readouts
(`integration/readouts/cn_*/arm.txt`, closed-arm / open-arm occupancy of every node's stage-B maps against the deposited closed
map and the consensus), the leaf scorecard (`integration/final/lsc/arm.txt`), site r_s5's ranking files, and a per-particle
profile of the root's four groups and of r4's two (`scripts/hetero_tree_pose_study/r4_profile.py`, output
`pose_study/r4_profile.txt` on batch2; one 4-core job of 1 min 45 s).

**Origin.** The root was seeded from site r_s3 (the strongest kept pole) plus two route-A extra states: site r_s4's pole became
stage-B state 3 (group 3, r3) and site r_s5's pole became state 4, which the choose routing delivered as group 4 = r4 (100,330
particles, the stage-B state 3 of the final numbering). r_s5 is candidate 3f2: a central site (centre 14 A from the molecule's
middle, 22 % of the molecule, 29,383 voxels above 0.5), and its kept pole is the merge of three of the four latent tails
(z1hi, z3lo, z3hi; the fourth, z1lo, became the rest). Three tails of two different latent directions sharing one deviation map
inside a central site is the shape of an outlier or quality axis, not of a conformation. All four root groups passed the re-test
(weakest pair z 27.7), the locality guard (1.01 to 1.06) and the junk screen; r4's 14 leaves are all `leaf`, none junk.

**Maps.** Every map in the r4 lineage has a weak arm in both positions: closed-arm occupancy 0.08 to 0.22 and open-arm 0.41 to
0.72 at the nodes r4, r41, r42, r411, r412, r421, r422 and all 14 leaves, called "intermediate" by the readout because both
occupancies are low, not because density sits between the positions. For comparison the corner r1311_s2 reads 1.00 / 0.58, r2's
open states 0.42 / 1.2 to 1.3, the root's own state 3 (r4's seed) 0.17 / 0.57. r4 holds 693 of the 8,720 closed-reference rows,
a share of 0.0069 against the base rate 0.0126.

**Particles** (root refinement, per group):

| group | n | turned > 30 deg | > 150 deg | median angle | median shift change | stage-B corr p10 / p50 / p90 | view TV vs rest |
|---|---|---|---|---|---|---|---|
| r1 | 252,458 | 0.065 | 0.027 | 3.7 | 0.72 px | 0.382 / 0.398 / 0.415 | 0.087 |
| r2 | 192,760 | 0.087 | 0.037 | 4.3 | 0.83 | 0.381 / 0.395 / 0.412 | 0.056 |
| r3 | 146,844 | 0.099 | 0.041 | 4.5 | 0.90 | 0.380 / 0.393 / 0.410 | 0.089 |
| r4 | 100,330 | 0.172 | 0.070 | 5.7 | 1.27 | 0.379 / 0.389 / 0.405 | 0.154 |

Defocus (1.28 to 1.30 um) and the odd/even balance do not differ. Turning is a low-SNR phenomenon in every group: of the
particles in the node's lowest stage-B-corr quartile 26 to 33 % turned, of the highest 0.7 to 1.4 %; the turned particles'
median corr is 0.381 against 0.396 for the rest, and 41 % of them turned by more than 150 degrees (near-flips), 27 % by 90 to
150. r4 differs from its siblings in two ways: its corr distribution sits a quartile lower (median 0.389, the node's first
quartile is 0.386), and within every corr quartile a larger share of its particles turned (0.33 / 0.11 / 0.04 / 0.014 against
0.27 / 0.05 / 0.02 / 0.007). Inside r4 the same pattern repeats one level down: r42 (36k) turns 6.2 % against r41's 2.2 % at
equal corr, with a slightly higher defocus (1.31 against 1.27). The closed lineage r1 turns 0.6 to 1.0 % per group.

**Verdict.** r4 is best read as a quality-and-alignment state, not junk and not a conformation: a group of lower-SNR particles,
gathered through a central-site pole that merged three latent tails, whose poses are the least stable in the tree and whose
maps lose the arm in both positions because the arm is what blurs first when poses scatter. It is not junk by the screen or by
defocus, and its corr deficit is modest; it is not a pose artefact in the micelle sense either, since the turned particles are
spread over all four groups and r4 merely holds more of them. It costs the tree 14 leaves and 100k particles (14.5 % of the
data), and it is a real-data instance of the audit's quality-axis failure mode (item 13) and of route A seeding from a site
whose pole is not a conformation. Two cheap checks would settle it, neither run here: reconstruct one r4 leaf from its top-corr
half and its bottom-corr half and read the arm (if the arm appears in the top half, the "intermediate" maps are a quality
blur); and see whether r_s5's pole survives the subspace site test, whose distinctness test should reject a pole whose
endpoints are three collapsed tails (the joint gate replay of the 692k will show whether r4 reappears). Both the pose flag
(dflip 0.106) and a stage-B-corr quality flag (r4's median at the node's first quartile) would have marked the root split that
created it.

**Half-split check (run; job 54855816, 16 cores x 2 min 40 s = 0.7 core-h).** Leaf r4111_s1 (9,061 particles; its leaf map reads
0.216 / 0.642, "open") was split at its median stage-B corr and each half reconstructed at the leaf's own stage-B poses (selection
on the merge-up's master project, which carries every leaf's poses; box_crop 88, no post-processing), both read through the saved
readouts' chain (`make_view.py --box 160`, a fixed raised-cosine low-pass, so the readout does not depend on a half's FSC):

| map | particles | stage-B corr p50 | closed-arm | open-arm | open : closed |
|---|---|---|---|---|---|
| top-corr half | 4,531 | 0.406 | 0.263 | 0.714 | 2.7 |
| bottom-corr half | 4,530 | 0.390 | 0.185 | 0.606 | 3.3 |
| full leaf | 9,061 | 0.398 | 0.216 | 0.642 | 3.0 |

The arm density follows particle quality in both positions: the top half, with half the particles, exceeds the full leaf in
both occupancies and the bottom half falls below it, while the open-to-closed ratio stays near 3 in all three maps. Neither half
resolves into a state (r2's open leaves read 1.2 to 1.3 open-arm, the corner 1.00 closed-arm), so the quality split exposes a
density gradient, not a hidden conformational subpopulation. The "intermediate" maps of the r4 lineage are a quality blur of an
open-arm population, which confirms the verdict above. The remaining check, whether r_s5's pole survives the subspace site test's
distinctness test, waits on the joint 692k replay. Scripts `scripts/hetero_tree_pose_study/r4111_halves.sbatch` and
`r4111_halves_flags.py`; data `pose_study/r4111_halves2/` on batch2.
