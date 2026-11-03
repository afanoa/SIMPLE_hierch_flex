# The 692k RNase P tree on the combined code (2026-11-02)

Run `integration/f692_lf_s1` on batch2 (job 54864239, 13:30 to 18:44 EST), build `math_reworks/build_lf` = branch
`claude/tree-lever-floor` at c5176bddc: the step-1 cleanup, the subspace site test, the SVD site discovery, the draw-free
single-pass split test, the took-off pole rule and the effect floor, on top of recipe17d's keys (`integration/f692_lf.sbatch` =
`f692_r17d.sbatch` with the build swapped: maxdepth 4, route A to depth 2, minpop 3,000, junk screen off, nparts 12, nthr 8,
nodepar 12, maps contrast, no rho floor). The reference is the recipe17d run `f692_r17d_s1` of 2026-10-29 (5 h 31 min, 873 core-h,
one closed-corner leaf r1311_s2). Readouts: `integration/readouts_lf/node_*.txt` (every node's closed-state content, read by a
follower job as each split landed), `integration/final_lf/` (leaf scorecard, embedding landscape), `route9/landscape/lf/`
(docking landscape against the 55 deposited apo conformers), memory `integration-692k-lf-launch-2026-11-02`.

## 1. Root and depth 1

| | recipe17d root | combined code root |
|---|---|---|
| groups | 4: 252k, 193k, 147k, 100k (r4, the quality blur) | 5: 163k, 146k, 136k, 128k, 120k |
| seeding site | candidate 7 (7f1) | candidate 7 |
| extra states (route A) | the two central sites 3f1, 3f2 | the peripheral sites 1 and 2 (poles kept at z 14 to 15) |
| split statistic | one-draw test | draw-free 0.915 (z 37.9); guard 1.00 to 1.07 on all ten pairs; floor P / N 248 to 361 |
| quality flag, max pair dflip | (r4: dflip 0.106) | off; 0.031 |
| closed-reference recall of the closed group | 0.56 (purity 0.019) | 0.38 (purity 0.020), 60 % of the closed rows in groups 2, 4, 5 |
| root wall, cost | 115 min | 107 min, 91 core-h |

The new discovery folds the two central sites into one whose poles are too weak to seed, so the r4 quality blur is not seeded at
the root; the closed particles spread over five groups instead of four. At depth 1 group 1 handed 89 % of its closed particles to
one child, r13 (46,991 particles, purity 0.063, closed arm 0.88), which recipe17d reached only at depth 3 (r131, 0.87).

## 2. The closed state (user verdict 2026-11-02: "yes they are all closed")

| leaf | particles | closed-reference rows (of 8,720) | purity | closed / open arm | on the deposited swing axis (matched) |
|---|---|---|---|---|---|
| r1312_g1_2 (r1, r13, r131) | 12,427 | 954 | 0.077 | 1.04 / 0.63 | +11.5 A, past every deposited closed conformer |
| r2112_g1_2 (r2, r21, r211) | 13,202 | 1,259 | 0.095 | 0.89 / 0.66 | +6.6 (unmatched axis) |
| r1321_g1_2 (r1, r13, r132) | 12,044 | 739 | 0.061 | 0.89 / 0.61 | +8.0 (unmatched) |
| r1311_g1_2, r1322_g1_2 | 12,470, 10,050 | 748, 520 | 0.060, 0.052 | 0.77 / 0.57, 0.71 / 0.45 | |
| r5111_g1_2, r5121_g1_2 (r5, r51) | 9,567, 8,159 | 415, 396 | 0.043, 0.049 | 0.68 / 0.80, 0.62 / 0.65 | |
| recipe17d corner r1311_s2 | 12,010 | 816 | 0.068 | 1.00 / 0.58 | +10.0 |

The scorecard's strict corner rule (closed >= 0.95 and open <= 0.6) admits none of the new leaves: r1312 misses it on the open arm
by 0.03. The user judged the 10 A unsharpened maps of r1312, r2112 and r1321 beside the old corner and called all three closed, so
the closed state is retained, in three lineages, with 2,952 closed-reference rows in closed leaves against 816 in the one corner
before; the best three-leaf union reaches F1 0.137 against the closed reference. The rule undercounts and must not be used alone
as the hard gate; the closed-arm occupancy and the closed-reference rows per leaf are the readout.

## 3. The rest of the tree

r3 (136k) is the pure open state throughout (closed arm <= 0.2 in all its leaves). r4's subtree r42 to r45 (94k, eight leaves)
has the weak-both-arms maps of the old r4 (closed 0.1 to 0.4, open 0.2 to 0.4); only r41 (34.5k, 938 closed rows) carries closed
content. No statistic flagged it this time (pose flip share 0.105 like its siblings, stage-B corr within 0.08 sd), consistent with
the pose-conditioned re-test of the same day: it is a population and quality difference, not an alignment artefact.
108 nodes, 66 leaves (recipe17d 84 and 67); the merge-up fuses 66 leaves into 63 units (recipe17d 67 into 63).

## 4. Landscapes

**Embedding** (`final_lf`, both trees' leaves on the new root's flex embedding): transport distance between the two leaf
landscapes 0.016 against 0.021 to 0.025 to the shuffled null; the new tree covers 93 % of the old one's extremes, the old 62 % of
the new one's; at the depth-4 cut 60 states against 46.

**Docking against the 55 deposited conformers** (`route9/landscape/lf`, the 2026-10-29 recipe; figures
`lf_landscape_armspan_dock_2026_11_02.png`, `lf_landscape_dock_2026_11_02.png`, `lf_landscape_dockmatch_2026_11_02.png`; in the
figures the script's old labels apply: "integration r17d" is the combined code, "old full692" is recipe17d):

| | deposited | combined code (66) | recipe17d (67) |
|---|---|---|---|
| swing axis s1, plain docking | -10.1 to +7.5 A | -11.1 to +9.4 | -10.9 to +8.2 |
| s1, amplitude-matched | -15.2 to +9.8 | -16.4 to +11.5 | -15.8 to +10.0 |
| 2 A bins of the deposited swing with a leaf | | 9 of 9 | 9 of 9 |
| best arm correlation to a conformer, median / min | nn 0.917 | 0.800 / 0.705 | 0.789 / 0.635 |
| conformers whose best leaf is in this tree | | 32 | 23 |
| conf43 (most open) best match | | r3121_s1 0.827 | r3111_s2 0.805 |
| leaves outside the deposited hull (joint MDS, plain / matched) | | 11 / 15 | 8 / 10 |

The entire swing is produced in both trees. The closed end is further in the new tree (r1312 at +11.5 A matched, past every
deposited closed conformer and the old corner's +10.0). The new tree's open extreme (-15.3 to -16.4 A) is the r4 quality subtree,
weak maps with near-zero open-arm occupancy, so the credible open edge is the same in both trees. Every deposited conformer has a
leaf of either tree within 3 A in the swing plane; no leaf is a strict copy of a deposited conformer in either tree (the pipeline
offset of 2026-10-29). The pairwise docked-correlation coverage is in section 6.

## 5. Wall and cost

| | recipe17d | combined code |
|---|---|---|
| wall | 5 h 31 min | 5 h 14 min |
| core-hours (step jobs under the run) | 873 (10,130 jobs) | 617 (12,805 jobs) |
| root phase | 1.91 h | 1.79 h |
| body | to 4.20 h | to 3.91 h |
| end phase | 1.32 h (merge-up 73 min) | 1.30 h (merge-up 78 min) |
| site layer core-hours | 332 (38 %) | 88 (14 %) |
| stage A + B core-hours | 279 | 249 |
| disk after pruning | 94 GB | 83 GB |

The 29 % saving in core-hours is the subspace site test: the per-site site_rank jobs (169 core-h in recipe17d) are gone. The wall
moved by 17 minutes because the root and the merge-up are serial; the merge-up port of the same day (branch
`claude/optimistic-chandrasekhar`, 65 min on the saved tree) and the core budget (about 20 min at the root) are the remaining
scheduling levers, as the audit's section 9 found.

## 6. Pairwise docked correlation and coverage

Every pair of the 188 maps (55 deposited, 66 new leaves, 67 old) rigid-docked at 10 A, amplitude-matched to the deposited mean,
box 64 (`docked_cc.py`; figures `lf_landscape_docked_cc_2026_11_02.png`, `lf_landscape_coverage_vs_count_2026_11_02.png`):

| | deposited | combined code | recipe17d |
|---|---|---|---|
| within-set nearest-neighbour CC, median | 0.978 | 0.942 | 0.952 |
| best CC to a deposited conformer, median / min | 0.978 / 0.953 | 0.922 / 0.825 | 0.912 / 0.746 |
| projected into the deposited-only MDS: hull area (deposited 0.208) | | 0.136 (65 %), 2 of 66 outside | 0.145 (70 %), 1 of 67 outside |
| joint MDS: leaves off the deposited cloud | | 48 of 66 | 49 of 67 |

The new leaves are closer to the deposited set than the old ones (median 0.922 against 0.912, and the worst leaf 0.825 against
0.746, since the old r4 mixtures are gone), cover about the same share of the deposited hull (65 against 70 %) and lie inside it
except r1321 and r2221, which sit just past its closed and open edges. The joint-MDS offset (48 of 66 off the cloud) is the
pipeline offset of 2026-10-29, a uniformly lower CC to every conformer, not new conformations. The worst-covered conformers are
again the large open ones (conf48, conf49: best leaf 0.928 against their own neighbour 0.995), so the coverage gap still tracks
atypicality, not the conformers' particle counts. The closed conformer conf18 has r1312 as its best leaf at 0.948.

## 7. Verdict

The combined code keeps the closed state (three closed leaves, the user's eye), produces the entire swing with the closed end
further out than the deposited set, matches the deposited conformers at least as well as recipe17d, loses the r4 quality blur at
the root but not its subtree, and costs 29 % fewer core-hours at the same wall. It is the production code line from here; the
merge-up port and the core budget are its next scheduling steps, and the quality population of r4's subtree the next open problem.
