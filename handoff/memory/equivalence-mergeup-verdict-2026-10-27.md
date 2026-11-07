---
name: equivalence-mergeup-verdict-2026-10-27
description: "Equivalence merge-up (merge only if shown equivalent + never fuse node-distinct siblings; test build ~/SIMPLE_smr2) on saved trees: keeps the small OLE states and full span, fuses 0 distinct siblings, reproduces MORE states across seeds; m33 closed pair kept apart (79 units); on nulls it keeps every false leaf (splits are the problem)"
metadata:
  type: project
---

2026-10-27, adaptive split/merge session (branch claude/adaptive-split-merge-b384e9, env SIMPLE_MERGE_RULE=equiv; outputs
ole_12707/split_merge_replay/):
- OLE 6_: delivered 14 units / 0 of 15 small leaves / 10 of 11 distinct siblings fused -> equiv 30 / 15 of 15 / 0.
  7_: 11 / 0 of 12 / 10 of 10 -> 24 / 11 of 12 / 0. base75, base90, MF, C1: 2x more units, 0 distinct siblings fused.
- Seed pairs (map counterpart >= 0.7 of noise ceiling): 6/7 reproduced units 10,8 -> 14,16; span along r_s2/r1_s4/r2_s4
  0.50/0.39/0.32 -> 0.55/0.52/0.34 (= leaves). base75 8,11 -> 9,10; base90 11,13 -> 20,18.
- cannot-link-only arms (cnl, cnl_z5) lose span and reproduce fewer units -> rejected.
- Truth: IgM m52M units 45 -> 67, merged units mixing truth-different leaves 22/26 -> 7/7; m69M 13/16 -> 1/2; qtoy2_bp 1/3 -> 0/1.
  Remaining bad merges come from pooled upper-level tests in diluted union masks -> fix = pairwise own-mask test (proposal R1).
- RNase P m33: delivered 32 units (closed pair r1111221_s1/_s2 fused) -> equiv 79, closed pair apart, 0 of 20 distinct
  siblings fused; base90 35 -> 69.
- qtoy null (amp_mstep_v1 qtoynull_base_s1, 5 false leaves): off 4 units, equiv 5. OLE null under
  SIMPLE_SPLIT_CONTRAST=yes: 1 leaf (root T 0.05 < 0.106; unset split at 0.81). qtoy null under yes still splits the
  root (T 0.09 vs 0.064) on no truth column (snrq AUC 0.50; unset 0.12).
- Caveat: "never overrule distinct" also keeps FALSE node splits; the split test has ~1-in-3 residual null false positive even
  with SIMPLE_SPLIT_CONTRAST=yes ([[split-test-accepts-quality-split-on-null-2026-10-27]]) -> null trees decide.
Related: [[focus-three-threads-2026-10-27]], [[discrete-splits-capture-the-span-states-that-matter]].
- RNase P HARD GATE HOLDS (14:45): on m33 the equivalence merge-up keeps the closed pair r1111221_s1/_s2 as two states
  (79 units; today's merge-up gives 32 and fuses the pair). Null: it keeps all 5 qtoynull leaves -> needs the split-test fix.
- NATIVE PORT (2026-10-28): /scratch/cluster_scratch/atifao/equiv_port/equiv_port.diff (state_mergeup only, default off;
  recipe4 = recipe7 base file). Native inputs == offline mergeup_inputs.py on qtoy_recipe5, qtoy2_bp (19 vs today 16),
  CAHRA sub_a (8 vs today 5: today fuses BOTH node-distinct pairs), IgM m69M (53 vs 33), qtoy null r4 (2 vs 1). Gap: cannot-link
  covers stage-B leaves only; sibling child nodes ending K*=1 (r1_g1_2/r2_g1_2) are not cannot-linked.
- 2026-10-28: DEFAULT ON for recipe12 (equiv_port_default_on.diff; SIMPLE_MERGE_RULE=off = old rule). OLE recipe4: unset == offline (17 units), off == delivered (7). Copies of the equiv/qnull work: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/{equiv,from_scratch}.
