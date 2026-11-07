# Adaptive split/merge for hetero_tree: hand-over (2026-10-27)

Proposal: `/Users/atifao/SIMPLE/.claude/worktrees/adaptive-split-merge-b384e9/doc/implementation_notes/planned/hetero_tree_adaptive_split_merge.md`
Replay scripts and outputs: same worktree, `scripts/hetero_tree_math_experiments/split_merge_replay/` (`results/`).
Nothing is committed and no native code was changed. One cluster job ran: 54383490 (norm, 16 CPUs, 2 min) in
`/mnt/beegfs/elmlund/afan/ole_12707/split_merge_replay/`.

## Main finding

The merge-up fuses sibling leaves that their own node had already shown distinct in the protected 16-10 Å re-test.

| tree | distinct sibling pairs | fused | particles |
|---|---|---|---|
| OLE seed 1 | 11 | 10 | 36k |
| OLE seed 2 | 10 | 10 | 38k |
| RNase P m33 | 20 | 16 | 158k |

- The RNase P fused pairs include the closed pair `r1111221_s1`/`_s2` (mutual best matches of conf19 and conf18).
  - Node test: z 5.4, in its own 2.9 % mask.
  - Fused at merge-up level 7: one 35-unit test, 17.7 % union mask, M from 11,900 null statistics.
- All 13 scorecard-scored fused RNase P pairs have local change ≥ 0.034, against the 0.028 deposited nearest-neighbour median.
- OLE `r12111`: the node test and the merge-up test had identical inputs. One said distinct (T 0.32 vs M 0.235), the other
  merged (T 0.23 vs M 0.339). At this size the test is a coin toss.
- OLE map check (map_compare metric, 15 Å):
  - The leaf sets reach the same extremes in both seeds. The merge-up units lose 6-25 % of the span along `r_s2` / `r1_s4` /
    `r2_s4`.
  - Reproduced states: 15 → 11 (seed 1) and 17 → 9 (seed 2).
  - The specific sibling differences mostly do not reproduce across seeds (median 0.27 of the ceiling). Interior partitions
    are stochastic.
- The 7 OLE K* = 1 pairs are "unresolved, not same", but splitting them is not supported: 7 of 10 pairs have the same
  counterpart in the other seed.
- m33's node splits still had the `em_min` 1,000 floor. 13 of its 23 K* = 1 nodes are drains with no re-test.

## Proposal

- Three verdicts (distinct / equivalent / unresolved). Merge only on equivalence at the consistency margin (the tree's smallest
  demonstrated effect). Never overrule a node's own distinct verdict.
- Pairwise merge-up tests in pair masks with a fitted threshold.
- On the saved trees this gives units = leaves on OLE and keeps the closed pair.
- `minpop` replacement (power to detect the parent's effect + Rosenthal-Henderson guard): specified, not testable offline;
  needs a native A/B.
