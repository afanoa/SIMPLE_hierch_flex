---
name: klem-six-conformer-700k-and-junk-flag-2026-10-19
description: "realistic benchmark: 700k halves, six conformers 46:42:49:27:47:29 (~7.4k per conformer per half) -> 7 states after the cross-half final pass (29+27 merged at 0.35, 42 split in two, one mixed); dtree_v11 = root-child junk flag (free-component share > 0.8 in both halves)"
metadata:
  type: project
---

2026-10-19. 700k synthetic halves (sim55_700k_A/B, refined poses lat_55_700k_A/B, run dirs synth_55_700k_A/B prepared by
synth_tree_prep_full with CHAIN=subset_only.sbatch). pick_confs.py chose 46:42:49:27:47:29 (four reference clades,
pairwise 0.34-0.59). dtree_v10 TAG=six: root K* 5, r.1 -> 46 / 42 / 42, other children leaves; cross-half final
pass (tree_finalize2, final_xem): 7 states = 46, 49, 47 clean-ish (0.68-0.73), 42 twice (0.59-0.73), 29+27 merged
(0.35 apart, unresolved although 46/42 at 0.34 split), one mixed 49/47/46 state. NMI conformer 0.31 -> 0.36.

dtree_v11.sbatch (= v10 + junk flag; v10 untouched): after a split at depth <= JUNKDEPTH (0), junk_children.py marks a
child junk when the node board put more than JUNKF (0.8) of its particles on the free components in BOTH halves;
junk children get junk.flag and are not split; tree_leaves.py appends " junk" in .names. Dry run on finished roots:
two_50k r.2, null_50k r.2, 55_200k r.3 at 0.99-1.00 (truth junk 0.95-0.97); the junk-free 55_700k root put 0.37/0.47
of real child r.1 on the free components, hence 0.8 in both halves rather than a 0.5 mean. End-to-end test: TAG=jflag
on two_50k (job 54008196). Related: [[klem-cross-half-final-pass-2026-10-19]], [[igm-like-synthetic-armsigm1-2026-10-19]].

CLOSED-STATE TESTS (dtree_v12 + cross-half final pass, same 700k halves): six + closed 18 at ~4k per half (8 %, sixc8):
5 states, the closed conformer gets its own state at 38-47 % purity (rest 42/46/47), but 47 merges with 49.
Six + 18 at 1000 per half (2 %, sixc2, MINGRP 600): 6 states, NO closed state (18 spread, <= 5 % of any child).
two_50k (49 open vs 18 closed, ~6k each + junk) separates them cleanly (0.92-0.98). Real data: closed ~2 % but
~7k particles per half; earlier it surfaced only by zooming a family-1 leaf ([[rnasep-closed-state-not-pure-in-klem-leaves]]).
