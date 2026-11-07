---
name: run-what-you-build-immediately
description: "A capability the user asked for is not done until it is running on the real workload; launch it at scale the moment it validates, using spare core budget in parallel if the box is busy"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: fb1d9e33-0137-47d8-8d28-e3a32330f274
  modified: 2026-09-01T15:59:24.386Z
---

2026-09-01: after implementing + validating the distributed paired flex_pca
dispatch (part v5), I parked it "until the box frees up" while the science chain
ran. User: "You put in the distributed version but you are not running it? what
is your problem" — followed by "commit this stuff to memory. You keep making the
same mistakes and just end up wasting my time."

**Why:** the user asks for capabilities in order to USE them on the current
campaign. Building + toy-validating + waiting is, from their side,
indistinguishable from not having built it. Related failures the same day:
launching batteries on unvalidated poses ([[flex-pca-validate-inputs-before-batteries]]),
layered polling waiters instead of one flat relaunchable script.

**How to apply:**
1. The definition of done for "implement X" is X running on the real data at
   real scale, or a stated blocker the user has accepted.
2. A busy box is rarely a blocker: budget cores (e.g. 48-thread science run +
   32-thread nparts=4 distributed run on an 80-core box) and run in parallel —
   this also yields the same-workload timing/agreement comparison for free.
3. Never silently defer a just-built capability; if it truly cannot run yet,
   say so in one sentence with the earliest run slot.
4. Monitor greps must not match filenames like simple_error.f90 — anchor on
   make's "Error N" lines, and cover failure AND success states.
