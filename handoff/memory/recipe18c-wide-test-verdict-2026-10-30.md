---
name: recipe18c-wide-test-verdict-2026-10-30
description: "recipe18c wide test (all new switches default off; all-off == recipe17d exactly): SUPPORT=order fails every RNase P closed gate (cuts the flexible arm), qtoy2 and the continuum; INHERIT keeps r111, loses r131's corner by .03, r1211 slightly worse; POSEGATE inert off micelle; ROUTE neutral-to-worse except discrete toys"
metadata:
  node_type: memory
  type: project
  originSessionId: 9643ce3d-fc7a-4dc1-860a-59db10867e7c
  modified: 2026-10-30T02:49:34.128Z
---

2026-10-29/30, integration session (recipe18c = recipe17d + inherit_on_r17d + route_v5 + micelle stack [tree_mskfile pair, fixposes,
SUPPORT, POSEGATE], hetero_tree_20261028/recipe18c; recipe18d = + guard fix locality_guard_children_r17d.diff d6e45d9f).
- All off == recipe17d line for line (L20_6k md1, identical labels).
- SUPPORT=order: closed arm 0.34-0.56 in every leaf on r131, closed_r111 and r1211 (FAIL); qtoy2 NMI .684 -> .532, closed F1 .95 -> .47;
  continuum eta2 .185 -> .027. Nulls pass (re-test fails). Only for micelle-type data, never general.
- INHERIT (+POSEGATE): r111 identical (2 corner leaves); r131 best closed leaf .92/.52 vs .97/.54 (corner lost) but F1 .147 vs .123;
  r1211 .642 -> .612 (open arm .61 -> .73); continuum ends reached. r131 seed 2 (recipe18d): OFF corner .99/.51, INH none
  (.92/.51) -> INH FAILS the corner 2/2 seeds (pole-less r1/r2: fflex seeds deliver the corner, inherited seeds do not).
- POSEGATE: root flip shares 0.04-0.33 % on RNase P, nulls, rare: never fires; on micelle it equals fixposes.
- ROUTE: RNase P neutral/worse (r1211 .642 -> .583-.596), continuum/rare/OLE worse; helps only discrete toys (mm3q 1/8 -> 6/8 combos).
- Guard skips: only where choose routes to EM groups (mm3q 3/14, qtoy2 INH 3/12, SUP rare runs); none in RNase P/OLE/IgM trees.

**Why:** the user (2026-10-29): "these finished threads need to be tested more widely"; "we need a method that works on everything".
**How to apply:** none of the switches is a general default; a combined method must pass RNase P closed gates first. Related:
[[continuum-inherit-verdict-2026-10-29]], [[micelle-fixture-pose-scatter-verdict-2026-10-29]], [[final-routing-pass-verdict-2026-10-29]].
- FINAL (2026-10-30, ~690 core-h, recipe18c_switch_matrix_2026_10_29.md): SUPPORT fail; INHERIT fail (r131 corner 2/2); POSEGATE
  harmless; ROUTE discrete toys only; guard fix (locality_guard_children_r17d.diff) judged all 5 former skips, all kept, partitions
  unchanged = ready to adopt alone. mm3q md2 fair baseline (recipe18d OFF) NMI .701: SUP level, INH slightly worse.
