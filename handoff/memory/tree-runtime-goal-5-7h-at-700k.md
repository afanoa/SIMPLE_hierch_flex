---
name: tree-runtime-goal-5-7h-at-700k
description: "STANDING (user 2026-10-25): a hetero_tree on the ~700k RNase P set must finish in 5-7 h wall so they can test / develop; the design features got 'really expensive for not sure how much benefit' -> optimize before adding features"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-25T17:22:01.456Z
---

2026-10-25 about 13:40, after I reported the new design's root at about 3 1/4 h (2x the baseline per node) and estimated about a day
for the whole 692k tree: "this has gotten really expensive for not sure how much benefit. Could we work on optimizing? Honestly I
think 5 - 7 hrs is the goal for 700k particle because otherwise I cant really test/develop".

**Why:** a tree that takes a day can't be iterated on. Every extra stage (every-site ranking, motion calls, the split-half gate,
null refits) has to earn its wall time, and the benefit evidence was still mixed (OLE reproducibility up, CAHRA worse from routing,
cdisc rare state lost to an older change).

**How to apply:**
- Measure wall time on the 692k tree as the primary cost, not just core-h. Report the critical path per node and where queueing
  adds to it.
- Before adding a feature, show its benefit on real data and its wall-time cost at 692k.
- Prefer subsample fits that score everyone once, node-size-scaled parallelism, and dropping redundant steps.
- Don't flood csbdevel with low-value test trees while a 692k run needs the cap.
- Related: [[dont-gate-tests-on-price-cores-idle]] (idle cores still apply, but wall time and benefit now gate design work),
  [[price-new-stages-lead-with-lean-scope]].
- 2026-10-28 PRECEDENCE: [[science-before-speed-2026-10-28]] — the 5-7 h goal applies to the tree that keeps the states; never trade a state-finding step for time.
