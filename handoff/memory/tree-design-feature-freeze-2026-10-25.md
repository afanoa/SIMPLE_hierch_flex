---
name: tree-design-feature-freeze-2026-10-25
description: "STANDING from 2026-10-25: no new hetero_tree design features until the usefulness scorecard (factorial session) shows benefit; speed work continues; the m33 root's chosen site 6 is mostly a per-particle QUALITY axis"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T18:01:30.676Z
---

2026-10-25, user via the coordinator, after viewing the m33 root's site-6 maps: "this is pretty minor. not sure if its doing anything useful. I think we need a better way to assess if we are doing something useful and if something is worth doing."

- Feature freeze on the tree design. No new design feature or selection rule from this session until the scorecard shows benefit.
  - The factorial session builds the scorecard: the 22 Oct deposited-conformer benchmark, effect size per split, toy truth, seed reproducibility and wall time.
  - The coordinator relays any rule to build. The next expected one is a knob-free quality screen for motion candidates.
- Speed work continues: build 83 (concurrent site ranking, strongest candidate first) and the base80 timing run.
- Routing (SIMPLE_SPLIT_ROUTE best/band): let the battery finish and report it, but build it into nothing.
- Scorecard result 1 (rnasep_klem/scorecard/m33_root_m79/): the root's chosen site 6 is mostly per-particle quality.
  - corr falls monotonically across the bins at matched defocus (AUC 0.32);
  - the rigid core's CC drops 0.94 -> 0.88, and the arm moves only 2.1 A.
  - The pole-rule margin ("strongest") rewards it over the real arm axes: site 1 (closed 0.49 -> 0.28) and site 2 (6 A).
  - So design80 and design83 were cancelled.

**Why:** every design feature must earn its wall time and show benefit; the "strongest" margin can favour quality over conformation.
**How to apply:** don't add or tune design features unasked. A design tree is judged by the scorecard, not by K* or margins.
Related: [[m33-tree-wall-profile-2026-10-25]], [[package-tree-tests-2026-10-25]], [[split-routing-em-vs-stageb-2026-10-23]].
