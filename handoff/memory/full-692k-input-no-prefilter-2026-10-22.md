---
name: full-692k-input-no-prefilter-2026-10-22
description: "User 2026-10-22 ~14:10: 'why aren't we running on the full 692k particle set?'. Real-data runs start from ALL 692,392 apo particles (rnasep_klem/f692/subset_all.txt), never from an earlier run's subset (643k f692gen seeds, 536k junk-free); junk is handled inside the tree"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T17:54:28.247Z
---

The full RNase P apo set is 692,392 particles: `rnasep_klem/f692/subset_all.txt`, 1-based master-row indices into
`rnasep_coupled/full692k_2d/rnasep.simple`, which has 932,496 rows.

What earlier runs used instead:
- The f692st tree (rflex13) started from 643,140 rows via `ROOTROWS=rnasep_klem/f692gen/seed_states.txt`. That earlier run had already
  dropped 49,252 particles; why is unconfirmed.
- I had proposed 536k, which also removed the user's 107k junk-like leaves.

**Why:** filtering a new run's input with an old run's output is circular. It also contradicts the per-node junk design: the tree must
find junk itself, and flag junk states rather than having them removed upstream.

**How to apply:** start real-data trees from subset_all.txt. Mention any exclusion explicitly and get the user's OK first. Related:
[[junk-detected-per-node-not-only-root-2026-10-22]], [[f692st-final-deliverables-2026-10-21]].
