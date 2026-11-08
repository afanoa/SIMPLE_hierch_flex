---
name: coordinator-arms-own-wakeup-while-sessions-wait
description: "2026-11-08 user: 'why did you decide to stop': 4 session reports sat unread ~2 h (12:52-14:19 EST) after my context was compacted; the app holds side-session notes after a failed turn or a long run of report-only turns. Arm a Monitor heartbeat (queue poll, <= 30 min expiry, re-armed) whenever cluster jobs are in flight."
metadata:
  type: feedback
---

On 2026-11-08 my last update went out at 12:32 EST. The context was compacted at ~13:27, and no turn of mine ran until the user
asked at 14:37. Meanwhile rnasep-misfit, sparse-ppca2, site-chain and extremes-repro had reported. Nobody acted on them, and two
sessions sat idle on multi-hour jobs with nothing to wake them.

**Why:** the user expects autonomous progress ("keep integrating at testing autonomously"). A coordinator that only waits to be
woken stalls silently whenever the app stops delivering notes.

**How to apply:**
- While any cluster job is in flight, keep a Monitor armed in the coordinator: a read-only `squeue --me` poll that prints each job
  that leaves the queue. Re-arm it on every expiry; the expiry itself is the heartbeat.
- Every brief tells the session to arm its own waiter on its jobs.
- On waking, read all queued reports before replying.

Related: [[job-monitor-batch2]], [[relay-intermediate-results-from-background-schedulers]].
