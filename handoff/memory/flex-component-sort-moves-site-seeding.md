---
name: flex-component-sort-moves-site-seeding
description: "TRAP 2026-10-25: m8c75t-89t sort flex components; site_rank reads FFLEX coords (base80's 1,3 = m33's z1, z10; replay pair 1,5; match rows by position). The sort was first blamed for the RNase P closed-state loss, but it is nearly harmless (0.989 vs 0.992 at 3 stA0 passes); the cause is the 1-pass stA0 root, see sta0-root-passes-decide-rnasep-closed-state"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-26T00:31:45.645Z
---

2026-10-25, closed-state loss in full692_base80_s1 (root on m8c80t), traced read-only (scorecard/closed_trace.*.out).
- Each flex fit at a node logs its own sort: base80 node_r flex "z1..zK were components 5 1 10 7 9 8 4 2 3 6", fflex
  "... 1 4 10 6 3 8 2 5 7 9". The fits are identical to m33's; only the order changed. Eigenvalues are nearly flat (4062 .. 2640).
- site_rank reads node/fflex/flex_pca_coordinates.txt (simple_hetero_tree.f90 site_submit), so the FFLEX permutation is the one
  that matters. The flex-log permutation was first used by mistake (it gave the wrong replay pair 2,9).
- Verified by column correlation (first 20k rows, r = 1.000000 each): base80 z1..z10 = m33 z1, z4, z10, z6, z3, z8, z2, z5, z7, z9.
  So base80's site_dirs=1,3 = m33's z1 and z10; m33's z1, z3 = base80 site_dirs=1,5.
- The two files number particles differently (row 2 is particle 3 in m33 and particle 2 in base80), but the rows are the same
  particles in the same order. Match rows by position.
  - The cause: under node_proj=yes (the default since 16e52fbdd) the id column holds node-local positions; m33's holds project rows.
    node_r/rows.txt is identical in both trees, and stree_port/closed_split.py maps rows through it.
- Checks launched 2026-10-25 about 19:05:
  - closedrep/base80_root_sd15 (54300397): base80's root from recB0 with site_dirs=1,5, on m8c90t;
  - base90 (54300456): a full 692k tree on m8c90t, base80's keys, nodepar 16.
- Seeding outcome:
  - m33 kept pole {z1lo, z3hi} (130k seeds, 18.7 %) against the rest.
  - base80 rejected pole 1 (60 % share) and kept z3hi (242k, = m33's z10) against the rest.
- Root K=2 split (re-test z 34.6 vs 3.71, PASS, so the 4e re-test did nothing here):
  - m33 kept all 11,231 rows of its closed leaves (r1111221_s1 / _s2, conf18 / 19) on one side, ending in a pure r1111221;
  - base80 split them 56 / 44, and they never exceed 6 % of any later node or leaf.
- base80's closed-arm leaf r1211111_g1_2 (0.84) has a sharp core but a mixed arm, so it fails the paired test: 0/3 closed recalled.
**Why:** the user's hard gate is "we need to retain the closed state"; a seeding axis change silently removes it.
**How to apply:**
- To replay m33's seeding on a sorted tree's files, use the fflex permutation (base80: site_dirs=1,5). A rerun that reads a
  sorted run's fflex output from disk keeps the sort even on an unsorted build.
- Check any build that touches flex component order this way. Score closed recall first on every tree (score_tree.sh, qtoy).
  Related: [[tree-scorecard-and-blur-quality-screen-2026-10-25]], [[rnasep-closed-state-must-be-retained]].
- DECISION (coordinator, 2026-10-25 about 19:00): the sort is REVERTED for tree builds. m8c90t = m8c89t without the sort (and
  without the Gamma pairing patch) is the build from here. Builds m8c75t-89t carry the sort, so their RNase P closed-state results
  are confounded (qtoy (b) / (c) included; (b') / (c') on m8c90t are the clean arms).
- CORRECTION (port, 2026-10-25 about 19:40): the sort is NOT the whole cause. sd15 (base80's recB0 + m33's exact directions,
  site_dirs=1,5) still splits the closed rows 6,782 / 4,449 (0.604); root ARI vs m33 0.023; pole z3hi only (223k); z1lo fails
  its row test (-2.76 SE).
  - Root inputs compared: mask.mrc and fflex states are identical in m33 / base80 / sd15 / base90. stA0 differs (m33's states
    5-8 overlap 0.79-0.90; base80 = base90 0.24-0.32).
  - The only difference in the stA0 command: maxits 3 (m33's ROOT, m8c33) vs 1 (77d50f3c5, 23 Oct, "ITA0 = 1", the user's speed
    decision from cdisc / CAHRA C11 toys; inside published d3513ff70). m33's body (m8c36) already ran 1 pass, as do qtoy (a)
    and (b').
  - Downstream, recB0's site table changes: m33 has 6 sites with top region "8"; base80 / sd15 / base90 have 7 sites, top
    region "1", with "8" 4th. The root's poles are ranked in site_01's mask, so they come from another region.
  - Being confirmed: an stA0 rerun at maxits=3, m33's root from its recB0 on m8c90t (54302084), and the 4th cell (3 passes +
    the sort's directions 1,10). Next test build m8c92t = 3 stA0 passes at the ROOT only, as m33 ran; base92 is the closed-state
    plus timing test. The default stays ITA0 = 1 until the user decides.
  - CONFIRMED 20:07-20:40:
    - base90's stA0 command at maxits=3 on m8c90t (closedrep/base90_sta0_its3) reproduces m33's root stA0: final-pass
      populations within 0.3 %, overlaps within 0.01;
    - base90's root (unsorted, one pass) equals sd15 to the last digit: z3hi-only pole of 224,824 seeds; closed rows 6,782 /
      4,449 (0.604); root ARI 0.023.
  - base92 = 54303755 on m8c92t: the snapshot-only patch .claude/patches/sta0_root_three_passes_on_m8c91t.patch.
  - 2x2 at the root (stA0 passes x seed directions; closed rows' largest share, root ARI vs m33):
    - 3 passes: (1,3) 0.992 / 0.953 (m33 replay on m8c90t: m33's exact pole, 129,678 seeds); (1,10) 0.989 / 0.935;
    - 1 pass: (1,3) 0.604 / 0.023 (sd15, base90); (1,10) 0.564 / 0.020 (base80).
    - VERDICT: the pass count decides the root; the sort's directions are nearly harmless once stA0 converges.
- LESSON: a speed change validated on toys (ITA0 = 1: same seeds and states on cdisc / CAHRA) changed the RNase P root. Gate
  speed changes on a cheap RNase P root replay (closed rows by child) too.
- The cause and the 2x2 are summarised in [[sta0-root-passes-decide-rnasep-closed-state]]. Whether the sort returns depends on the
  OLE sort-off pair (16_ / 17_ on m8c90t vs the sorted base75 12_ / 13_, both 1-pass roots).
