---
name: sta0-root-passes-decide-rnasep-closed-state
description: "CAUSE 2026-10-25: the RNase P closed state is lost at the ROOT when stA0 runs ONE refinement pass (ITA0 = 1, 77d50f3c5); with 3 passes (m33's root) it stays together (0.99) whatever the seed directions. 2x2 verified; base92 (3 root passes only) is the whole-tree test; default is the user's call"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-26T00:31:26.121Z
---

2026-10-25 evening, port replays scored on m33's closed rows (r1111221_s1 + _s2, 11,231 rows) at the 692k root.

| stA0 passes | seed dirs | run | closed rows by child | share | root ARI vs m33 |
|---|---|---|---|---|---|
| 3 | (1,3) | m33 itself (root on m8c33) | 11,231 / 0 | 1.000 | 1.000 |
| 3 | (1,3) | m33 recB0 replay, m8c90t | 11,146 / 85 | 0.992 | 0.953 |
| 3 | (1,10) = the sort's | m33 recB0 replay, m8c90t | 11,113 / 118 | 0.989 | 0.935 |
| 1 | (1,3) | sd15 / base90 | 6,782 / 4,449 | 0.604 | 0.023 |
| 1 | (1,10) | base80 | 6,339 / 4,892 | 0.564 | 0.020 |

- 77d50f3c5 (2026-10-23, "stA0 runs one refinement pass (ITA0 = 1)") was the user's speed decision. It was based on fast-tree cdisc /
  CAHRA C11 toys, where 1 pass gave the same seed site, seeds and stage-B states.
  - m33's ROOT ran 3 passes (m8c33); its body (m8c36) ran 1 and still isolated the closed state at depth 7.
  - qtoy (a) on m8c36 runs 1 pass at the root too.
- At the root, 1 pass is not converged:
  - an stA0 rerun at maxits=3 on m8c90t matches m33's stA0 within 0.3 % / 0.01;
  - after pass 1 it sits at base90's states, and passes 2-3 move them (state 8: 36k -> 81k).
  - The unconverged stA0 changes recB0's site table: m33 has 6 sites with top region "8"; 1 pass gives 7 sites with top
    region "1". The root's poles are ranked in site_01's mask, so they are ranked in the wrong region.
- The flex component sort is nearly harmless once stA0 is converged ([[flex-component-sort-moves-site-seeding]]).
- Test: m8c92t = m8c91t + .claude/patches/sta0_root_three_passes_on_m8c91t.patch (ITA0_ROOT = 3 at depth 0 only, snapshot-only).
  base92 (54303755, base90's keys, nodepar 16) runs it, with base90 (1 pass, otherwise identical) as the control.
**Why:** the user's hard gate is retaining the closed state ([[rnasep-closed-state-must-be-retained]]). A speed change validated on
toys silently broke it on real data.
**How to apply:**
- Don't run RNase P trees with a 1-pass root.
- Gate every speed change on a cheap RNase P root replay (closed rows by child, root ARI vs m33) as well as on toys.
- Whether 3 root passes become the default is the user's call. A convergence stop for stA0 (repeat until the populations stop
  moving) is the knob-free alternative, raised for after the design freeze.
