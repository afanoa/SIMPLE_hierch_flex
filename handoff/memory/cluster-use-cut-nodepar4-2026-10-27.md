---
name: cluster-use-cut-nodepar4-2026-10-27
description: "STANDING (user via coordinator, 2026-10-27 ~22:00): cut cluster use; ask the coordinator before launching anything new; every run nodepar <= 4 (supersedes the nodepar 16+ default)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 01619764-7985-49aa-b2e9-539d9cd5a16f
  modified: 2026-10-28T01:38:54.821Z
---

User instruction relayed by the OLE/adaptive-partitioning coordinator session on 2026-10-27 late evening: "we don't need that
many cores". Caps are unchanged (csbdevel 642, norm 600 per user), but use far less.
**Why:** the user saw the account pinned at both caps with ~3,900 CPUs pending across sessions, most of it low-value queued work.
**How to apply:** do not launch new trees or batteries without asking the coordinator first; any tree launched uses nodepar <= 4
(this supersedes the nodepar 16+ default in [[tree-launch-fast-defaults-2026-10-25]]); tree controllers request 2 CPUs / ~16G;
hold seed-2 replicates until seed 1 shows a signal. Related: [[efficient-testing-over-wide-queues-2026-10-27]],
[[keep-cluster-busy-queue-next-job-first]] (now overridden: do not pre-queue).
