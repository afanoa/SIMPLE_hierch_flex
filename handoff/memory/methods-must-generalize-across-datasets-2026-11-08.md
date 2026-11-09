---
name: methods-must-generalize-across-datasets-2026-11-08
description: "USER 2026-11-08: 'this method has to be generalizable and have an advantage across data sets' -> every tree change is judged on a panel (truth fixtures + RNase P + OLE + EMPIAR-10076), never on RNase P alone; the near-tie seeding failure already generalizes (RNase P r1, OLE depth 1)"
metadata:
  type: feedback
---

2026-11-08 ~23:30, after the robustness redesign was proposed (guide 7.10 thread 14): "and of course this method has to be
generalizable and have an advantage across data sets."

**Why:** a fix tuned to RNase P is a knob, not a method; the paper needs an advantage that holds on several molecules.

**How to apply:**
- Evaluate every tree change on a panel:
  - truth fixtures, where diversity is countable: qtoy2, qtoyflex, mm3q, the IgM 10 % / 5 % fix and nulls;
  - RNase P (rare closed state; 55 deposited conformers);
  - OLE / EMPIAR-12707 (seed reproducibility on record);
  - EMPIAR-10076 (deposited assembly classes EMD-8440..8449).
- Metrics: reproducibility across seeds, the number of distinct honest conformations (diversity), recall of deposited states,
  and nulls staying one leaf.
- The motivating failure is general. Near-tie seeding decisions diverge runs on OLE at depth 1
  ([[ole-repro-seeding-near-ties-2026-10-28]]) and at RNase P r1 (guide thread 13), so fixes target the mechanism, not a molecule.
Related: [[replicate-or-redesign-positives-only-2026-11-08]], [[flex-pca-must-be-knob-free]],
[[discrete-splits-capture-the-span-states-that-matter]], [[landscape-must-be-reproducible-and-splits-data-derived]].
