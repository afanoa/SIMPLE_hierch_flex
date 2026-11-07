---
name: klem-null-leaves-sort-by-alignment-quality-2026-10-19
description: null-set leaves (one conformer) are told apart by alignment quality (score/pose error), not view coverage; view matching leaves 14-15 states; fix under test = per-shell spectral match before the difference FSC
metadata:
  type: project
---

Null 200k set (conformer 49 + 20 % junk), tree + consolidation with VIEWMATCH=1 (leaves thinned to the whole set's
100-bin direction histogram, vm_merge_test.sbatch): the data-derived merge still leaves 15 -> 14 states (floor 12 A)
and 15 (auto band scan); remaining difference-band FSC median 0.54 / 0.86. View coverage was NOT the separator.
leaf_diag.py (per leaf vs random subsets): every leaf is junk-free (0.00 of 20 %) and well aligned (angerr 4.1-4.8 deg
vs 11.4 whole set); leaves differ in refinement score (corr 0.42 vs 0.45, z +-30) and pose error (4.1 vs 4.8 deg).
A leaf of worse-aligned particles is a BLURRED copy of the same map; blur differences are reproducible between even
and odd halves, so the difference-map FSC calls them real. Known-K (18+49, 50k): 8 states, purity 0.98, no merge.
Fix under test: SPECMATCH=1 in merge_data_derived.py (every leaf's full map scaled per shell to the consensus signal
spectrum via the half-map cross spectrum, then the deviation/difference FSC); specmatch_test.sbatch.
Also: squeue on batch2 needs -u atifao (not afan); 4M stack merge now streams (merge55.py, 48 GB).

**Why:** the merge must be invariant to per-state resolution/blur or quality-sorted duplicates survive on every real set.
**How to apply:** run null + known-K after any merge-rule change ([[klem-scaling-synthetic-55-2026-10-18]]);
view matching alone is not sufficient ([[klem-region-local-pooling-2026-10-17]]).

Update 2026-10-19 noon: the in-board `specmatch 1` (patch_specmatch.py) only reaches the single-process solve in
simple_commanders_klem_board.f90; the DISTRIBUTED master (simple_klem_distr.f90 klem_board_master, own Wiener solve ~l.413)
ignores it, so the nparts=10 A/B gave a byte-identical latent. Region-board cells of two_50k_A showed no corr/pose-error
split (cells_nuisance, but its unweighted spread is dominated by tiny cells: low power). Leaf level (two_50k): pose error
4.2-5.2 deg and corr 0.41-0.43 across leaves, defocus flat. Merge-level SPECMATCH ported to merge_cross_tree.py
(cross_tree_test.sbatch SPECMATCH=1 XS=_sm), runs 54002793/4.
