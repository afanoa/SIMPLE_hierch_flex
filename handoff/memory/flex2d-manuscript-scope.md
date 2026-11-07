---
name: flex2d-manuscript-scope
description: flex2D paper is flex2D-only with the FULL theory kept; cls_split lives in a separate draft (decided 2026-10-03)
metadata:
  node_type: memory
  type: project
  originSessionId: 95245f2d-90ec-4676-978e-28a5533c3f54
  modified: 2026-10-04T03:37:34.266Z
---

2026-10-03 length cut: the user wants the flex2D manuscript (`doc/manuscripts/flex2D_cls_split_2026_10_03/main.tex`,
IUCr format) to be about flex2D only, and to KEEP the full flex2D theory (no tightening of Section 2). cls_split goes
in its own draft, `doc/manuscripts/cls_split_2026_10_03/main.tex` (theory, sibling test, tab:cs-sets, per-data-set
observations; intro and discussion are placeholders). The main paper mentions subclassing only as the way l_0.5
is measured for the label-band rule.

**Why:** the user found the combined paper too long. cls_split adds no state information beyond the view, so it
weakened the flex2D story.

**How to apply:** do not re-add cls_split results, figures or panels to the flex2D paper. Montages with a cls_split
panel are kept under `figures/with_cls/`. Supplementary material (Tables S1–S2, Figures S1–S5) is one-column after
`\onecolumn`. Pre-cut backup: `main_v8_prelengthcut.tex`. Related: [[judge-2d-partitions-by-heterogeneity]],
[[flex2d-findings-2026-09-21]].
