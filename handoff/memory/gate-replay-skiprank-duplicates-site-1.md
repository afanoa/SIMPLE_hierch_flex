---
name: gate-replay-skiprank-duplicates-site-1
description: "TRAP: stree_gate_replay2 SKIPRANK=1 on a source tree that ranked only site 1 makes the last site a second copy of site 1 (fsite_<node>_s<last>), so a replay's 'two passing candidates' are one motion"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-25T10:15:32.407Z
---

Found 2026-10-25 in the m33 design replay (stree_port/gate/m33_design, job 54228924).

- SKIPRANK=1 writes the site count into `sites/seed_site.txt`. The controller finds no record in that site's dir, so `start` goes to `advance(S_RECB0)`, which always ranks `site_01.mrc`. It writes into `seed_dir`, which is `fsite_<node>_s<last>`.
- Result: candidates 1 and <last> have byte-identical `motion/c*/bins.txt` and identical motion stats. Only the gate score of <last> is its own, because the gate refit reads the site mask.
- In real trees `seed_site.txt` changes only through `next_site`, which ranks the matching mask. A scan of 4,239 nsite/fsite dirs (2026-10-25) found every `pcg_mskfile` matching its site number. The restart path would mislabel only a k > 1 site step that lost its record. A guard was proposed to the coordinator, not applied.

**Why:** a replay read showed "the strongest of 2 passing" with both candidates the same site.
**How to apply:** in a SKIPRANK replay, count only sites whose rankings the source tree really made. Check the `md5sum` of `motion/c*/bins.txt` before calling two candidates distinct.
Related: [[package-tree-tests-2026-10-25]].
