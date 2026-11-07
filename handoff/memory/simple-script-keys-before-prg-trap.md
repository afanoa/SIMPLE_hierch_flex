---
name: simple-script-keys-before-prg-trap
description: "TRAP: a cmdline key set before prg= is written first into the job script and the job dies 'prg= or test= must come first' (exit 1); hit twice in hetero_tree patches (crop keys 2026-10-27, rowmap in fs_rank 2026-10-28)"
metadata:
  node_type: memory
  type: feedback
  originSessionId: c8c6f555-6cba-44a4-9ec4-48197b959f2f
  modified: 2026-10-29T01:21:40.498Z
---

SIMPLE job scripts write cmdline keys in insertion order. Any `cl%set(...)` before `cl%set('prg', ...)` produces
`simple_exec rowmap=... prg=site_rank ...`. That dies at once with "ERROR! prg= or test= must come first on the command line"
(exit 1) while the controller carries on.
- 2026-10-27: crop keys in deliver_state_maps.
- 2026-10-28: rowmap in recipe15's fs_rank. Every site-layer rank failed silently; the node path still looked right, because those
  ranks are optional inputs to state_seeds.
The fast gate does not exercise either path.

**Why:** the failure looks like a missing optional result, not a crash, so equivalence checks can pass by accident.
**How to apply:** in any new cmdline code set prg first. After a patch, grep a fresh step's hetero_step.sh for the key order, and
check every site/step folder's hetero_step.exit is 0 before comparing results.
