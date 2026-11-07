---
name: locality-guard-skips-em-group-splits-2026-10-29
description: "TRAP: state_mergeup locality_guard keeps a split UNTESTED when choose routes to EM groups and every stage-B state's majority falls in the same EM group (log shows only 'K* 2 -> 2', no L/L_ref line)"
metadata:
  node_type: memory
  type: project
  originSessionId: f6f31e73-5d9d-4609-b06e-14cf14e50c04
  modified: 2026-10-30T00:38:52.816Z
---

Found 2026-10-29 on qtoy2_inh_md2_s1 node r4 (holo, 6,721 rows). Re-test z 3.61 PASS; choose -> EM groups (5,074 / 1,647).
Both stage-B states' majorities landed in EM group 1, so smap(s) == smap(t) for every state pair. locality_guard then computes no
pair ratio, rmin stays huge, and K* is kept. The result was two 92 %-holo leaves plus a 76 %-junk sliver. Present in recipe17d
(simple_state_mergeup.f90, locality_guard pairing loop); not fixed.

**Why:** the guard is the only check against diffuse, quality-like splits since the rho floor was dropped
([[floor-off-guard-verdict-2026-10-29]]); a silent skip lets exactly those through.
**How to apply:** when reading a split log, a "locality guard: K* a -> b" line with no "L / L_ref" line before it means the guard did
not test the split. Count these before trusting guard verdicts. A fix needs a group-level map comparison when no state pair spans
the groups. Related: [[tree-continuum-fixture-and-verdict-2026-10-29]], [[null-split-diffuse-real-split-localized-2026-10-28]].

FIXED 2026-10-29 (route9/guardfix/locality_guard_children_r17d.diff, md5 d6e45d9f, against recipe17d / 5709bc083; gated):
- Unspanned group pairs, and missing stage-B halves, are judged on the delivered children's own half maps.
  state_split reconstructs them and pads them to box 88; same pair_locality, same 0.5 bar.
- Children-path calibration: nulls 0.22-0.36, real splits 0.87-1.03.
- The known qtoy2_inh r4 instance reads 0.639 and is kept.
- The 692k r17d tree was never affected.
