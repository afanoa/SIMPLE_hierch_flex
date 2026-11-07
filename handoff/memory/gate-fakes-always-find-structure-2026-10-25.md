---
name: gate-fakes-always-find-structure-2026-10-25
description: "OLE design pair 10_/11_: every fake-direction refit's posteriors match a real ranked start at |rho| 0.03-0.99, far above 3/sqrt(N) (0.01-0.016), so a 3-SE 'found structure' drop leaves no fake null (only the permutation null); seed 2's r11 promotion miss is the halves rule (every group member >= 3 SE both directions), not a bug"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T19:20:49.487Z
---

2026-10-25, for the coordinator, from CAHRA's private_split_check of ole_12707/10_tree_design_s1 vs 11_tree_design_s2 (m8c77t).

Floor:
- Stored motion/null/d*/scores.txt rho (column 3) = max |Spearman| of a fake refit's posterior to the node's real starts.
- Kept fakes: s1 0.029-0.891, s2 0.238-0.885; dropped ones >= 0.9. All 18 kept fakes exceed 3/sqrt(N) (0.0106 at r, 80,599 rows; 0.014-0.016 at r1/r2).
- So site_rank's EM finds the node's real motions from any start, and the fake null is not noise. The floors were s1 0.269, s2 0.742.
- A knob-free 3-SE drop leaves only the permutation null (4 random-fifth gate recs per floor node), which neither tree computed.

Promotion miss (seed 2 r11):
- The c2 axis is in site 1's posteriors (z1lo / z1hi, +0.825 / -0.809 to seed 1's c2), but merge_poles_halves accepted only z3hi.
- z1lo failed at -1.07 / -0.90 SE; z1hi passed the row test at 3.14 but not the column at 0.95.
- The single pole held every seed (share 1.0), so pole_rule rejected it; "0 0" = rejected.
- r111 site 1 failed by a hair: z3lo's column test at 2.97 and z1lo's row at 2.77 sank both groups.
- Ran on m8c77t after the resume (fix A in place).

Follow-ups, 2026-10-25:
- Floor: ADOPTED as the permutation null in m8c87t. OLE floors are 0.025 / 0.020 (stree_port/floorrep), against the fakes' 0.269 / 0.742.
- Pooled group test, offline (stree_port/pooled_halves3.py: an equal-norm proxy from the saved craw + quarter matrices).
  - At s2 r111 it keeps both c2 ends (3.14/3.97, 4.71/5.13 SE; +0.858 / -0.681 to s1 r11 c2). At r11 it does nothing (no merged groups).
  - Extra groups: OLE 0 / 5, 1 / 3, 1 / 2 per pair; m33_s1 17, base80 8, pkg_s1 10.
  - 11 of 12 OLE extras reproduce at |rho| >= 0.66 in the other seed's LINEAGE (node, parent, grandparent). Same-node matching
    understates this, because diverged seeds' same-named nodes hold different particles.
  - Proposed as an off-by-default key computing the exact pooled maps; awaiting the coordinator.

**Why:** both results bound what the gate and the site rule can promise.
**How to apply:** don't call fakes a noise null. Treat "candidate missing" as the halves rule's strictness before suspecting bugs. Any change waits for the scorecard ([[tree-design-feature-freeze-2026-10-25]]).
Related: [[package-tree-tests-2026-10-25]], [[ole-seed-split-choice-2026-10-25]].
