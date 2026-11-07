---
name: klem-cross-half-final-pass-2026-10-19
description: "final pass with a CROSS-HALF EM (splitbench4.py --xem, tree_finalize2.sbatch) beats the own-half soft EM on every tree: two_50k NMI 0.44-0.55 -> 0.63-0.64; 55_200k 13 -> 10 states, clade NMI 0.21 -> 0.25-0.26; a drained ~300-particle group sets the re-test's equal N"
metadata:
  type: project
---

2026-10-19 (user: "try the final pass to cross-half re-assignment"). toy55/splitbench4.py = splitbench3 + --xem N:
after the first decision, each half is scored against the OTHER half's group maps (consensus + deviation x cross-half
FSC per shell, <= 20 A, 480 pixels/particle, soft M-step, prior from own weights), groups under --em_min dropped, then
the same re-test. tree_finalize2.sbatch writes dtree_<MODE>_<TAG>/final_xem${XS} and final_nmi.py compares leaves /
old final / new final. splitbench3.py stays unchanged (running dtrees use it).

| tree | states old -> new | NMI conformer (or clade) leaves / old / new |
| two_50k | 5 -> 5 | 0.53,0.42 / 0.55,0.44 / 0.64,0.63 (converges in 4 iterations) |
| null_50k | 4 -> 4 | 0 (one conformer); real in junk states 0.011 -> 0.007 |
| 55_200k lat_or | 4 -> 4 | clade 0.53,0.50 / 0.53,0.50 / 0.54,0.56 (17 % still flipping per iteration) |
| 55_200k lat | 13 -> 10 (5 real + 5 junk) | clade 0.23,0.22 / 0.21,0.21 / 0.25,0.26 |

**Why:** the own-half soft EM is frozen by self-inclusion ([[klem-tree-representation-assessment-2026-10-19]]).

**How to apply:** use --xem 6 in final passes. TRAP: a group drained to ~300 particles became the re-test's equal N
(276), cutting its power: 13 -> 10 was mostly that. With --em_min 1000 (final_xem_min1000): 55_200k lat 12 states
(6 real + 6 junk; old 6 + 7), clade NMI 0.25/0.27, conformer NMI 0.10; two_50k 4 states (2 real + 2 junk), NMI 0.64.
A posterior-sum drain check still left a 367-particle MAP group, so the check now uses MAP counts and
tree_finalize2 defaults EMMIN=1000 (not yet rerun).
