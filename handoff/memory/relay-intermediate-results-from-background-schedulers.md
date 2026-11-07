---
name: relay-intermediate-results-from-background-schedulers
description: "a background scheduler that prints tables mid-run does not notify me until it exits: the 2-pole and CAHRA tables sat unread for 1.5 h (2026-10-23); make each table its own background task, or check the output on a timer"
metadata:
  node_type: memory
  type: feedback
  originSessionId: a4258345-b249-48ff-a0d3-d33e99ae33ac
  modified: 2026-10-23T17:00:36.497Z
---

2026-10-23, fast-tree settings study. My Mac-side scheduler (tset_sched.sh, one long background task) submitted 25 jobs and
printed each set's truth table as the set finished. The task notifies me only when it EXITS, so the 2-pole and CAHRA tables
(11:26 / 11:29) were never relayed. The coordinator had to ask at 12:58. The scheduler also stalled before scoring the last set.

**Why:** people downstream were waiting on those tables (the user's morning decision on the refinement schedule), and a silent
background task looks like progress while nothing reaches them.

**How to apply:** when a background job produces results in stages, make each result its own background task (it exits and
notifies when the table is ready), or schedule a check of its output. Never run a multi-hour, multi-result loop as one task and
wait for its exit. Related: [[keep-cluster-busy-queue-next-job-first]].
