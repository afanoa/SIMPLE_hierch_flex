---
name: tree-data-derived-lp-verdict-2026-10-29
description: "difference-SNR band: safe at 11 roots but FAILS deep: r131 subtree loses both child splits and the closed-corner leaf (2026-10-29); stays off"
metadata:
  type: project
---

Replay 2026-10-29 (hetero_tree_20261028/dsnr_band, dsnr_54604768.jsonl), on saved recB/stageB half maps.
- RNase P closed splits are NOT 10-15 A: the m33 closed pair has 20-12 A SNR_d 0.25 and I90 21 A; the user's one-helix sub-pair
  (recipe16 r111 r1) is the finest closed pair (20-12 SNR_d 1.25, I90 16.7 A).
- IgM arms (I90 21 A) are not coarser than RNase P; OLE is the coarsest (24.5 A).
- Nulls put 32 % of their difference in the split re-test band (14.7-10.2 A); real splits 3-19 %.
- Native replay (recipe17b + SIMPLE_TREE_LP=dsnr, 11 roots, 305 core-h, seeds identical): no gate lost; bands 10.6-21.3 A.
  The qtoy null is refused by the split test itself instead of by the guard. Closed r1211 F1 .528 -> .547, qbp arm unchanged.
  Other roots unchanged. Off by default (feature freeze). Diff in the shared notes folder, tree_dsnr_band_on_r17b.diff.

- DEEP NODE (r131 closed subtree of the 692k recipe17d integration tree, fixture at its depth-3 context: kway_maxdepth=-1,
  split_ntests=460, maxdepth 1): off keeps the closed corner (leaf r2_s1 0.97/0.54, 4 leaves). dsnr keeps the root split (rho .091)
  but both children fail z at stage-B 18-19.5 A (2.30, 1.32 vs 3.71): 2 leaves, NO corner leaf. The hard band is NOT safe at deep
  nodes: the split test's 14.6-10.2 A re-test needs the classifier to use that band. Ancestral shrink (shelved by the user): the root
  is refused by the 0.03 rho floor alone (z 14.5 passes, rho .0298), 1 leaf, no corner. 82.6 core-h.

**Why:** the user's premise ("heterogeneity at different resolutions") holds only weakly; the size coarsening remains.
**How to apply:** use I90 from dsnr_band.py before claiming where a split lives. Related: [[tree-data-derived-lp-verdict-2026-10-27]].
- DEEP-NODE TEST 2026-10-29 (r131 closed subtree, 44,334 rows, recipe17d + dsnr, 0.03 floor as launched; ancestral_dsnr/build_ad):
  root kept (z 20.84, rho .091, guard 1.03) but children r1/r2 run stage B at 18.0/19.5 A and fail the re-test (z 2.30/1.32 vs 3.71;
  off arm z 6.07/5.25 at 10 A) -> 2 leaves, NO closed-corner leaf (off arm: r2_s1 0.97/0.54). Cause: the band removes the 14.6-10.2 A
  re-test band from classification. Same coupling sank ancestral shrink at the r131 root. Hard band = not safe deep; keep it off.
