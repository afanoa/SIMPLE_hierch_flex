---
name: full-runs-only-for-meaningful-change-2026-11-08
description: "USER 2026-11-08: 'I dont want a safety check ... I dont want to run it if its not going to meaningfully change results' -> no full 692k trees as baselines or safety checks; launch one only when a lever has shown a clear gain on node replays and is expected to change the RNase P result"
metadata:
  type: feedback
---

2026-11-08 ~17:00. I had proposed, and the user had approved, a 692k wip-defaults tree that I expected to come out "same as
production". When I said so, the user cancelled it: "I dont want a safety check. We had like 8 closed leaves in the presentation at
least. also I dont want to run it if its not going to meaningfully change results."

**Why:** full trees cost ~600-900 core-h and a night. The user wants runs that can move the RNase P closed-state result, not
confirmations that nothing broke.

**How to apply:**
- Before proposing a full 692k run, say what result it is expected to change and the evidence for that: a lever with a clear
  gain on node replays (r1, r13, r131) with seeds and honest reads.
- Don't propose 692k runs as baselines, safety checks or default flips.
- Count production's closed leaves the way the user does: closed-arm >= .6 plus the user's eye. lf has 7, not the strict 3.
Related: [[wip692-full-tree-vs-production-2026-11-08]], [[test-on-subtrees-not-full-692k-2026-10-28]],
[[do-the-requested-job-not-diagnostic-sprawl]], [[integration-lf-closed-three-leaves-2026-11-02]].
