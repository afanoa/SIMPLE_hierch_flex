---
name: klem-region-local-pooling-2026-10-17
description: region-local merge of the 250 zoomed f692c sub-states, consolidation at cuts 0.4/0.5/0.55 (32/55/67 pooled states, 5.7-8 A), closed pair survives only from cut 0.55; the cut is NOT data-derived
metadata:
  type: project
---

Region-local merge (merge_states.py REGMASK=1: deviations multiplied by the leaf's own region mask, found at
fam?/regions{,_motion}/region_N.mrc) of the 250 FSC-filtered zoomed sub-states of f692c: 55/79/128/191/244 states at
similarity 0.5-0.9 (the whole-map tree gave 3/5/15/63/210: consensus-dominated). Consolidation (consolidate.sbatch,
TREE=merge_tree_zoom_reg.npz LEAFGLOB=zoomregion* CTAG=_zoomreg MINP=1000; one-job wrapper pooled_cut.sbatch +
pooled_analysis.sbatch) reconstructs every group from all its particles: cut 0.5 = 55 states 5.7-11 A (mostly 6-8),
cut 0.55 = 67, cut 0.4 = 32. The two pure closed sub-states (occupancy 0.92, 4.8k ptcls) are each other's NN at 0.74 and
merge with four intermediates (0.75) at cut 0.5; from 0.55 they stay a state of their own (state_24, 7.6 A). Strict
recovery of the 55 deposited: 4 (cut 0.5), 7 (cut 0.55), 12 (unpooled zoom). User asked "is this 55 data derived": no,
fixed threshold, coincidence with the paper's 55. Proposed data-derived rule: merge only while the pooled-map difference
does not reproduce across halves (same test as the leaves). Delivered: ~/Downloads/rnasep_klem_20261015/
deposited_conformers/pooled_states_55 and pooled_states_67_closed_kept (maps, README table, galleries).
Traps: consolidate.sbatch's inline gallery inherits LEAFGLOB (0 maps) - run pooled_analysis.sbatch instead;
consolidate_states.py falls back to zoomregion dirs for zoom names; region-name regex must accept region_1 and
regionmotion_1. See [[klem-leaf-states-are-redundant-across-regions]], [[rnasep-closed-state-not-pure-in-klem-leaves]].

**Data-derived merge (merge_data_derived.py, late 2026-10-17):** rule = two states are one unless the FSC of their
difference half-maps (scale-matched deviations, lpc cancels) is >= 0.143 averaged over the 1.25 x floor .. floor band,
inside a soft mask; greedy, least-resolved difference first; labels npz -> consolidate_states.py. Lessons: first-crossing
resolution is useless for difference maps; whitening not needed; masks must be low-passed with a REAL cutoff (Gm: 5 % at
15 A; the merge_states G is 60 % at 8 A and leaks); own FSC of deviation maps is inflated at high res by the common
-lpc term (ignore it). Sphere mask = whole-box test: the two closed zoom sub-states read "same" (0.11), 45 of 250 zoom
merges at 12 A, 122 at 10 A, 0 of the 72 leaves. Region-local soft masks (union of the two states' regions) = the
sensitive test: NO pair of the 72 leaves nor of the 250 zoomed sub-states is unresolved at 12 or 10 A (the closed pair
differs reproducibly inside the arm region, 0.69); the leaves are all distinct conformations at this depth. Also found:
node_states2 caps every leaf at --nmax 30000 (zoom 10000) most-confident particles (state mass up to 75k) -> UNCAP=1 in
consolidate_states.py takes every argmax-cell member; runs: f692c/f692e leaves + f692e zoom (jobs 53921133-5).
f692e zoom: 279 sub-states, closed 0.96 twice. Unrefined-quarter closed test: unrefined_closed.sbatch (53921090/1).
