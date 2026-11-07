---
name: dont-gate-tests-on-price-cores-idle
description: "STANDING (user 2026-10-25): 'not sure why you are held up by prices. we have 1200 avail cores for us and I assume they are mostly idle' - launch tests, report prices, don't hold them for approval"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-25T05:22:02.526Z
---

2026-10-25, about 01:45. I held whole-tree tests (OLE two seeds about 680 core-h, truth toys about 335-670, RNase P top levels about 430)
for a price go, and asked the port to cut the site-ranking cost first. The user: "not sure why you are held up by prices. we have 1200
avail cores for us and I assume they are mostly idle."

**Why:** the cluster capacity (about 1,200 cores for us) is mostly idle. Waiting on price approvals costs the user wall time and
attention, and spends nothing that's scarce.

**How to apply:**
- Launch the tests the work needs, in parallel, and report what they cost after (sacct). Don't ask for a price go first, and don't
  pare a test's method down to save core-h.
- Still keep the waste checks: price roughly to set sane -t limits for tests, stop runaway or failed jobs by id, avoid duplicate runs.
- This supersedes the 150 core-h-style approval lines in [[price-new-stages-lead-with-lean-scope]] for tests. Production runs keep
  their generous limits.
- Related: [[batch2-400-core-limit-and-worker-timeouts]], [[use-the-full-cpu-cap-for-tests]].
