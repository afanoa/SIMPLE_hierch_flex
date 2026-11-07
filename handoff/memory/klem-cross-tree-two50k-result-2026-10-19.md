---
name: klem-cross-tree-two50k-result-2026-10-19
description: "cross-tree (gold-standard) merge on the two-conformer 50k pair gave 8 states, not 2: every leaf is the well-aligned junk-free subset split by alignment quality (corr), a REAL reproducible partition, so independent trees cannot collapse it; fix must be nuisance-blind (in-board spectral equalisation / per-particle band)"
metadata:
  type: project
---

2026-10-19 10:30 (klem_xtree.54001346): two_50k halves (25k each, conf49 56 %, conf18 24 %, junk 20 %): 12 + 12
region leaves (all s01), all 12 matched (deviation corr 0.85-0.97, i.e. every leaf looks like the consensus), cross-tree
difference FSC median 0.27 at the start, only 4 merges -> 8 states; purity 0.7 (conf49 only), conf18 at 0.52-0.56 in
two small states. leaf_diag: every state has angerr 4.5-5 deg vs 12 for the set (z -5..-25), junk 0 (z -40), and the
states differ in image corr (z +11 vs -20): quality partitions. Same mechanism as the null set
([[klem-null-leaves-sort-by-alignment-quality-2026-10-19]]): a sharper cell map attracts well-aligned particles, a
blurrier one the badly aligned, in BOTH halves and BOTH trees, so the blur level is a latent axis and gold-standard
halves cannot remove it ([[klem-cross-tree-merge-design-2026-10-19]] does not work as a merge rule alone).
Candidate fixes (in the board, not post hoc): equalise every cell map's amplitude spectrum to the pooled one each
M-step (post-hoc SPECMATCH on leaves was too late); or a per-particle band/noise level so a blurry particle cannot
distinguish sharp from blurry cells. Test on the leaves first: leaf_bfac.py (relative B-factors, deviation
correlation before/after amplitude matching).

Null pair (synth_null_50k A/B, one conformer + 20 % junk, 2026-10-19 11:26): the same cross-tree merge gives 5 states per
half (12/11 leaves -> 5 at cut 12), purity 1.00 trivially. So the gate yields 5 false states at zero heterogeneity and 8 for
2 conformers: the count is dominated by the alignment-quality partition, not by conformers; specmatch A/B is the fix under test.
