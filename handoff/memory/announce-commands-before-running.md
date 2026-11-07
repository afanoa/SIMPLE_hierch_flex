---
name: announce-commands-before-running
description: "State in plain text what a cluster/long-running command will do BEFORE issuing it, and batch test runs into one submit-wait-report command instead of polling"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e5141833-c145-408a-ba48-7c3a3ffe0fbd
  modified: 2026-09-13T02:06:59.425Z
---

Say what a command is going to do before running it — especially SSH/cluster commands that
submit jobs, build, or write run directories. The user interrupted two tool calls in a row on
2026-09-12 for this; the second time the only complaint was that the command had not been
described.

Also: do not burn turns polling. Batch submit + wait + report + the comparison into ONE command,
and shrink the problem size for plumbing gates rather than running full production jobs.

**Why:** the user watches the tool stream and cannot tell an expensive or destructive command
from a cheap one until it has already run; repeated poll/check turns read as wasted time.

**How to apply:** before a cluster command, list in the response what it generates, submits,
writes, and compares, and note explicitly if it deletes anything. Then issue one combined
command. See [[flex-pca-ab-testing-traps]] and [[run-what-you-build-immediately]].
