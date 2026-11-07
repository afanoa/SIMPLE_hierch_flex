---
name: site-climb-lever-no-go-2026-11-02
description: NO-GO for the audit-6.5 rare-state lever (one hard E-step from the top-N* of the site EM posterior before the subspace test); 0 of 8 endpoints reach the bar at 6.5 % over 1-3 passes, and the pass destroys the pure 8-10 % endpoints and hands kept poles to the null roots; N* is an SNR optimum (3-10x the class), not a class size
metadata:
  type: project
---
2026-11-02, session tree-lever-floor (branch claude/tree-lever-floor, note hetero_tree_reworks_2026_11.md section 9.2; replays
math_reworks/lever/sites/c0-c3_* on build_lever, scored by lever/site_score.py against detlim/fix truth):
- 6.5 % site 5 (n_eff 116, bar r 0.275): endpoint r .13-.16 -> .11-.21 after 1-3 climb passes (best z 2.2); AUC unchanged .75-.79.
- 8-10 % site 1: the pass turns the rare endpoints (AUC .89, class 1.7-3.5k at purity .74-.79) into anti-poles (AUC .15-.17) and
  kept poles of 19-28k seeds at purity .03-.19 (the wrong end at 8 % s1 / s2 pass 3); 6.5 % nuisance site 3 gets a kept pole.
- Nulls: qtoynull / olenull root sites keep a 27k / 26k pole after one pass (none before) -> would seed null splits.
- Also rejected: the complement form of (B) (took-off class -> 1 - posterior): on the r131 seed-2 root site it kept the complement
  pole (4,566) and dropped the closed-population pole {z1lo,z3hi} (2,981 -> 349, rejected).
**Why:** the hard start at N* (15-40k) holds the rare state at purity .09-.12 vs the node's .065-.08, less contrast than the EM's
own class, and the prior N*/n (.4-.67) pushes the class past half the node; the posterior AUC is at the per-particle limit.
**How to apply:** do not propose rescoring/climb passes at the site as a rare-state lever; the lever would be per-particle
information (band, mask, pose). Related: [[rare-state-detection-limit-2026-10-29]], [[math-reworks-verdicts-2026-11-02]].
