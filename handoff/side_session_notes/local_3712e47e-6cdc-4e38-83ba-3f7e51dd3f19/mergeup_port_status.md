# Merge-up port (session tree-mergeup-port) — DONE 2026-11-02 16:30 EST

Branch `claude/optimistic-chandrasekhar-10506d` = `claude/tree-lever-floor-bccc29` (c5176bddc) + c569be829 (the port) +
55ee9c815 (the note, section 11). Nothing pushed; working tree clean. Fast-forward tree-lever-floor-bccc29 to 55ee9c815.

## Verdict
- Gates (math_reworks/gates/mu_* vs bc_*): the seven trees identical to bc in leaves, merge-up units (qtoy2 8 -> 8,
  r131 4 -> 4 both seeds, OLE 8 -> 8 with 0 merges), closed corner (r1_s2 1.00/0.49 seed 0; r2_s2 0.97/0.53 seed 1),
  NMI .812, OLE landscape W1 0.000 vs bc (coverage 1.00/1.00), nulls one leaf; merge-up walls the same (4 / 2 / 2 / 17
  min). Inside the merge-ups no pair reaches the equivalence test (weakest pairs z >= 6.1); bc's production test had
  reached it on r131 (seed 1: shown equivalent, kept only by the cannot-link; now distinct at z 12).
- 692k merge-up alone on the saved r17d tree (integration/mu_mergeup_r17d_s1, job 54866127, 48 cores, shadow
  split_stats): 64 min 44 s vs 73 min 4 s, 51.9 vs 58.5 core-h, MaxRSS 23.3 vs 22.4 GB; 65 units vs 63. Both runs merge
  r2311_g1_2+r2312_g1_2 and r2122_s2+r2511_s1; the saved run's other two merges (r3112_s2+r3111_s1, r2211_g1_2+r2212_g1_2)
  passed the old margin by 0.000006 / 0.000001 only at the pair's member-count n under the max-of-20-shuffles band
  (0.117); on the draw-free null (sd 0.008 vs 0.023, threshold 0.041) they are distinct at z >= 3.7. The four
  closed-lineage leaves are in no merge in either run. Whether those two pairs are one state each is a question for the
  maps.

## Optional like-for-like (not launched: 54864239 still queued, csbdevel 422 cores at 16:15)
After f692_lf_s1 finishes (its nodes carry split_stats.txt natively, no shadow):
  sbatch --export=ALL,TREE=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/integration/f692_lf_s1/1_hetero_tree,BUILD=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/build_mu/install,OUT=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/integration/mu_mergeup_lf_s1 -J mu_mergeup_lf /mnt/beegfs/elmlund/afan/hetero_tree_20261028/integration/mu_mergeup.sbatch
(48 cores, 30 GB, ~65 min; compares against that tree's own mergeup/merged_leaves.txt.)

Queue: no 692k step pending behind my jobs at any check (14:13, 15:25, 16:15). All my jobs are done.
