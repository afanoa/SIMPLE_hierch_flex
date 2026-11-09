---
name: production-tree-has-per-node-state-evidence-files
description: production f692_r17d_s1 tree: every node's stageB/prob_state_evidence.txt holds the refinement's own per-particle state posteriors (nats); softmax = responsibilities, sharper and more closed-informative than the hard labels
metadata:
  type: project
---

`$H/integration/f692_r17d_s1/1_hetero_tree/node_<id>/stageB/prob_state_evidence.txt` (r17d build; header "iptcl state T
evidence(1..nstates) in nats", -1e9 = no candidate; evidence = log sum over the neighbourhood candidates of exp(-T (d - dmin)),
T = pi wsqsum / (2 pftsz) x mask share). Read 2026-11-08 with node_levers/nl_resp.py (`self` labels): median max-responsibility
.86-.92, 43-54 % of particles > .9, a tenth near-uniform; closed rows sharper (r13: median .989); the closed-leaning child's
responsibility ranks the closed rows better than its hard label at every level (r .600 -> .631, r1 .682 -> .746, r13 .697 -> .804,
r131 .554 -> .580 AUROC). Maps at those nodes were pruned (only r1311 kept); b.simple projects remain.

**Why:** the user asked for responsibilities out of the refinement; they already exist for production at zero cost, while a
cross-build pass is invalid ([[wip-refine3d-redraws-poses-each-iteration-2026-11-08]]).

**How to apply:** for soft-label tests on the production tree start from these files; soft maps need a responsibility-table
importer into a state weight set (reconstruct3D m_estimator=flex), ~30 lines, not yet written.
