---
name: site-layer-off-for-test-runs-2026-10-27
description: "STANDING (user 2026-10-27): hetero_tree site_layer=no for development/test runs (A/Bs, null controls, timing trees); site_layer ON for deliverable OLE trees (per-site motion description)"
metadata:
  type: feedback
---

User, 2026-10-27, on the porting session's finding that OLE's factorial site layer adds ~1.4 h of sequential tail: "Turn
site_layer off for development and test runs (A/Bs, null controls, timing trees), but keep it on for deliverable OLE trees where
we want the per-site motion description."

**Why:** test runs judge leaves, speed or false splits, which the site layer doesn't change; deliverable OLE trees need its
per-site motion axes (the reproducible r_s2 / r1_s4 / r2_s4 motions used to judge landscape coverage).
**How to apply:** pass site_layer=no on every test/A/B/null/timing tree; pass site_layer=yes (or leave the default) only on trees
meant as OLE deliverables; say which one a run is when launching it. The 692k RNase P recipe already runs site_layer=no.
Related: [[efficient-testing-over-wide-queues-2026-10-27]], [[tree-launch-fast-defaults-2026-10-25]].
