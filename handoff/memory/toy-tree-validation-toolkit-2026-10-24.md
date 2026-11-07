---
name: toy-tree-validation-toolkit-2026-10-24
description: "Scorers for hetero_tree toy/CAHRA trees vs truth (routes/toy55, uncommitted in worktree laughing-carson-4e8ca1): validate_tree3.py (+units_prep.py for merged units), track_node.py, compare_spread.py, split_support.py, validate_conf.py; TRAP: validate_tree.py/2 counted the tree's 'junk' leaves as genuine"
metadata:
  node_type: memory
  type: reference
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-24T17:57:27.089Z
---

All in scripts/cls_expansion_proto/toy55 (synced to batch /mnt/beegfs/elmlund/afan/routes/toy55). Outputs go to
rnasep_klem/factorial/validation/.
- validate_tree3.py NAME TREE TRUTH sites|arms (validate_trees3.sbatch, env TREES 'name:kind:dir[:table]'), one SUMMARY line:
  - junk P / R, where removed = no leaf plus leaves named 'junk';
  - per-junk-group isolation and junk in genuine leaves;
  - purity, NMI, A-pure share, closed leaf, site-layer pole (factorial/table.txt);
  - every split's truth support (>= 0.2 share or 0.1 eta^2 = supported, weak, none).
- Merged units: units_prep.py turns hetero_tree_merged.txt (unit per row) into a flat layout and lists, per node, whether its
  leaves collapse to one unit. Then validate_tree3 scores the units: units.sbatch for sites, units2.sbatch for a kind per tree.
  E-step states map through the leaf -> unit map with --labels.
- track_node.py REF NODE TREE TRUTH: follows a reference node's rows down another tree to the first split that divides them, e.g.
  p27cd1's null split r211.
- compare_spread.py: tabulates SUMMARY lines, marking each arm within / beyond the p27cd1-vs-cd2 seed spread, or with --arms
  flagging BIG changes vs one baseline (IgM).
- split_support.py: split recurrence across seeds, see [[split-recurrence-across-seeds-2026-10-24]].
- validate_conf.py: the CAHRA (STAR conformation) version, open leaf and open-majority pool.
- TRAP: validate_tree.py / validate_tree2.py count leaves the tree names 'junk' as genuine leaves. The m8c17 IgM baseline's junk recall
  read 0.554 there; it is 0.917 with the junk leaf counted as removed. Use validate_tree3.py.
- Toy truths: simsites rnasep_klem/factorial/junk_test/simsites_truth_jt.txt (jtype 1 / 2, 20k each); IgM
  factorial/armsigm1_truth/arms_truth.txt; CAHRA EMPIAR_downloads/13991/data/source.star.
Related: [[m52-toy-arms-merged-units-2026-10-24]], [[tree-junk-pole-reentry-2026-10-24]].
