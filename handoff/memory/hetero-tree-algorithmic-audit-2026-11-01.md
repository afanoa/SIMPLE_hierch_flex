---
name: hetero-tree-algorithmic-audit-2026-11-01
description: "hetero_tree audit (user asked 2026-11-01): ~16.2k lines in the run path, 6.5x in 11 days; 26 % of the controller is the off-by-default motion layer; site_rank runs ~45 EM fits per site; split test builds maps twice and mixes null/real N; proposals in doc/implementation_notes/hetero_tree_algorithmic_audit_2026_11_01.md"
metadata:
  node_type: memory
  type: project
  originSessionId: 207d8490-49e6-4e5f-b348-619d653de945
  modified: 2026-11-02T00:10:29.904Z
---

2026-11-01: the user asked for an algorithmic audit of the tree ("extremely slow", "bloated", "site ranking / masking could be
more theoretically sound"). Written to /Users/atifao/SIMPLE/doc/implementation_notes/hetero_tree_algorithmic_audit_2026_11_01.md
(untracked, main checkout). Audited version: worktree funny-pike-0046ee (claude/tree-dist, uncommitted) on 21fac34db.

Key numbers: hetero_tree.f90 7,610 + state_mergeup 4,657 + landscape_sample 1,254 + commanders 1,057 + state_tree 709 +
tree_wiener 625 + UI ~250 = ~16,200 lines. Controller buckets: motion/gate/floor 1,842 (off by default), scheduler 1,511,
factorial site layer 1,191 (on by default, description only, 516 core-h at 692k), discovery+masks 939, step machine 821,
pole text-file bookkeeping 476. site_rank = 4 endpoints x 7 EM fits + 6 (rest) + 6 per kept pole (~45 per site).
Proposals: z-scored cross-half covariance for the mask core (no 0.2 x max); cross-validated SVD for site factors; one
subspace/principal-angle test on the cross-half deviation maps instead of the halves/quarters refits; map-space angle
instead of posterior Spearman for "same variable"; score-at-fixed-maps instead of EM climbs; a single split-test pass
(allm once, one FSC curve, one null at neq, numbers not log text for the merge-up).

**How to apply:** every proposal is a same-seed replay on qtoy2, r131, the nulls and OLE first; never a new default.
Related: [[tree-design-feature-freeze-2026-10-25]], [[rnasep-closed-state-must-be-retained]], [[m33-tree-wall-profile-2026-10-25]].

UPDATE 2026-11-01 (step 1 done, uncommitted in worktree funny-pike-0046ee on claude/tree-dist 9005a8c77, user builds): the
motion/gate/floor layer, SIMPLE_TREE_REENTER and the SIMPLE_SPLIT_{QSTRATA,EMNULL,KTGAP,RHOMIN,ORTHO,HPMASK,ROUTE,LOCALITY}
and SIMPLE_MERGE_{RULE,DELTA,DISTINCT} branches are gone (-2,260 lines; hetero_tree 5,725, state_mergeup 4,395). TRAP: builds
from this code ignore SIMPLE_SPLIT_RHOMIN / SIMPLE_SPLIT_ROUTE=em / SIMPLE_MERGE_RULE=off in old launchers; SIMPLE_SPLIT_CONTRAST=maps
is still honoured and still required by recipe17d. The math reworks run in session tree-math-reworks (worktree off claude/tree-dist).

UPDATE 2026-11-02 (user approved): claude/tree-dist (worktree funny-pike-0046ee) is at c7ac3b0bd = 9005a8c77 + step-1 cleanup
(e12be4e8d) + subspace site test variant D (26c350abb) + single-pass split test with numeric split_stats.txt and quality flag
(9c88430cd) + draw-free pair statistic (9ee5d5c72) + SVD site discovery as default (e567106fb) + variance-spread log line
(38cac2553) + notes, with state_split calling test_once (draw-free) as the DEFAULT (c7ac3b0bd). Not built, not gated jointly,
not pushed: the draw-free test and the SVD discovery were each gated alone; the combination needs its own replay (qtoy2, r131 both
seeds, three nulls, OLE) before any production run. Open after the night: endpoint scoring pass before the reproducibility test
(rare-state lever), guard + P/N >= 2 effect-size floor (oversplits), pose artefacts, merge-up on numeric stats + deletion of the old
test/decide path, site-layer climbs/sameness (6.4/6.5), step 2 opt-ins, the 59,145-seed anomaly in the new site_rank.
