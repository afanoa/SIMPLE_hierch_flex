---
name: do-the-requested-job-not-diagnostic-sprawl
description: "2026-10-04 the user stopped a chain of 7 control diagnostics with \"not sure why you are making a mess\" — run the requested job, report its one control result, ask before spawning investigation side-jobs"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T23:04:21.341Z
---

**What happened:** the user asked to align the 281k deselected 10180 particles to the 46k consensus. I instead launched a growing series of control experiments: greedy/prob/neigh × bands × grids, then a cross-dataset test on RNase P. The user interrupted: "I am not sure what the problem is. Can you align the deselected particles to the 46k consensus. Not sure why you are making a mess."

**Why:** the user wants the requested job done. A failing control is worth one sentence and one decision, not an open-ended investigation with new jobs and scripts each turn.

**How to apply:**
- Run the asked-for job in one script and report its single built-in check.
- If something looks wrong, report it once with one proposed fix and wait.
- Don't chain new diagnostic jobs without the user's go-ahead, even when they said "investigate". Keep that to 1–2 targeted tests, then report.

Related: [[small-fast-cluster-steps]], [[act-on-the-literal-ask]].
