---
name: null-acceptance-majority-leaf-2026-10-27
description: "STANDING (user 2026-10-27): null controls pass if the LARGEST LEAF holds >= 80 % of particles (one leaf not required; 50 % was a misread)"
metadata:
  type: feedback
---

User, 2026-10-27: "yeah idk if we are ever going to hit one leaf. But it would be adequate if the majority is kept in one leaf".

CORRECTED minutes later: "I meant more like 80%+ in one leaf lol. So its not good enough" -> bar = >= 80 %.

**Why:** a perfect single leaf on homogeneous data is not a realistic bar for the split test; the user wants the bulk of a
homogeneous population kept together.
**How to apply:** report the largest-leaf share (and merged-unit share) for every null tree; pass = >= 80 % in one leaf. Still report what
drives the minority splits (quality / view / featureless) and, for site_states, any extra site-pole states on the null.
Status at the time: qtoynull_base_s1 (d3513ff70) 53 % FAIL; qtoy null + contrast=yes FAIL (root 31.8k/28.2k); OLE null + SIMPLE_SPLIT_CONTRAST=yes = 1 leaf PASS.
Related: [[null-control-fixtures-2026-10-27]], [[split-test-accepts-quality-split-on-null-2026-10-27]].
