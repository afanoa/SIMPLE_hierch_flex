---
name: klem-leaf-states-are-redundant-across-regions
description: the chain's landscape count (count_states on focused cells) overcounts: a merge tree on the leaf state MAPS (cross-half deviation similarity, whole map) gives 5/8/14/32/52 states at 0.5-0.9 for f692c vs 60/83/110/184; regions of a family rediscover the same particle groups
metadata:
  type: project
---

2026-10-17: merge_states.py on the 73 level-1 leaf maps of f692c (even/odd reconstructions, scale-matched deviation from the
consensus inside the envelope, cross-half similarity, average linkage): nearest-state similarity 0.8-0.97 for most states,
self reproducibility median 0.99. Counts: 5 / 8 / 14 / 32 / 52 states at similarity 0.5 / 0.6 / 0.7 / 0.8 / 0.9, versus the
chain's cell-based count 60 / 83 / 110 / 184 at 0.6-0.9. The region boards of one family partition the same particles under
different focus masks, so the same groups reappear in every region, and the cell-based count compares deviations masked
differently per board (cross-board similarity ~0 -> inflated). The user had judged by eye that the leaf maps were blurry
near-duplicates needing consolidation; the map-based tree agrees.
Consolidation: consolidate.sbatch (cut -> pooled particles -> per-family reconstruct3D from all particles -> gallery),
outputs <run>/consolidated_<cut>/, galleries figs/gallery_<run>_consolidated_<cut>_pNN.png.

Consolidated f692c at cut 0.7 (13:17): 12 real merged states of 18k-87k particles, FSC 5.0-8.4 A (even/odd), 9 of the 10 docked conformers as best matches (46 missing), plus two residues of 272 and 4 particles with FSC 14.6 / 117 A (the data-support floor shows up by itself). f692d merge tree 6/9/17/33/53 (vs f692c 5/8/14/32/52): refinement does not change the consolidated count.

f692d (unrefined) consolidated at 0.7 (15:13): 11 real merged states >= 26k particles (+ 7.6k, 2k, 0.9k), FSC 4.7-9.4 A, 8 distinct docked best matches (46, 48 missing), 0 closed / 2 intermediate; f692c (refined): 12 real, FSC 5.0-8.4, 9 distinct (46 missing), 1 intermediate. Refinement's gain at the consolidated level is marginal (slightly better FSC and one more conformer).

USER VERDICT 17:30: the 73-86 leaf maps 'look pretty good... picking up most of the landscape', far more diverse than the 12 merged states; the merge over-collapses. Junk among the leaves = 12 of 73 with < 1000 reconstructed particles or FSC > 12 A, all with cross-half reproducibility 0.97-1.00 (junk reproduces). Fixes installed: NFREE_SUB=1 free component in child/feedback/region boards (junk sink) and a leaf support filter (own half-map FSC <= FSCMAX 12 A and >= MINREC 1000 particles) in merge_states.py / merge_by_particles.py. Particle-overlap (Jaccard) between leaf states of different regions is only 0.1-0.3: by particles they are NOT duplicates although their maps correlate 0.8-0.97 (consensus-dominated). Landscape = leaf states after the FSC filter, merged only on real particle overlap.

Soft state maps (soft_state_maps.py, 17:50): soft map = mean of the member cells' written cell_*.mrc (exact); half-maps = full + placed (D_half - D_mean) from klem_halves.npz, placement convention verified on a cell DIFFERENCE (corr 0.999): centred cube, phase origin at the box centre ((-1)^(h+k+l)), scale 1; D0/D1 are deviations from the board's OWN Wiener consensus (not toy_consensus.mrc), so never compare them against the consensus file. Per-particle weights = summed responsibilities over the member cells (weights_NN.txt) = stored weights for the M-estimated refinement (user: soft reconstructions default).

ZOOM ANALYSIS (19:26, f692c, 97 zooms -> 274 sub-states): 250 pass the half-map FSC/count filter; map-similarity tree 3/5/15/63/210 at 0.5-0.9; particle-overlap tree: NO overlap between any two sub-states (250 at every Jaccard cut); STRICT RECOVERY 12 of 55 (3 6 7 16 26 27 28 31 33 40 49 54) vs 8 for the level-1 leaves and 3 for the whole-map-merged 12. The zoom is the right direction; the landscape unit = zoomed sub-states filtered by FSC, merged only on particle overlap (none here).

USER 19:50 on the 12 conformer-matching zoomed maps: 'really bad resolution, once you go to so few particles you lose helix definition' and the closed one is not among them (its 0.92 sub-states match the closed GROUP without a 0.05 margin between closed conformers). Fix in progress: region-local merge of the zoomed sub-states (merge_states.py REGMASK=1: each leaf's deviation multiplied by its region's mask before the cross-half correlation; sub-states of one region that agree inside the region pool; union reconstructed) -> expected 30-60 states of 10-40k particles at 5-7 A with the purity kept.

**Why:** the honest landscape count is the map-based merge tree, not count_states; report it with the cut.
**How to apply:** replace step (5) of the chain by merge_states + consolidate; keep count_states only as a diagnostic.
