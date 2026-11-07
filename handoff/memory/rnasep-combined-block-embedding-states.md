---
name: rnasep-combined-block-embedding-states
description: "2026-09-20: the user's preferred RNase P readout = whole-molecule states clustered from ONE combined block embedding; BW_n1 (3 regions x 1 PC = 3 columns) and BW_n2 (2 PCs each = 6 columns) on 932k particles: all columns reproducible, 8 GMM states with region-localised, half-reproducible deviations (0.15-0.7 rel energy, e/o 0.86-0.96), fit ~35 min each; low rank per region is the right shape for clustering (30 columns drowned the placement)"
metadata:
  type: project
---

Runs `/mnt/beegfs/elmlund/afan/rnasep_regions/BW_n1`, `BW_n2` (launch_flex.sh with `pcg_mskfile=G10cov_labels_B3.mrc
neigs=1|2 rec_states=yes npreimages=8`, deconvolution on: 6 min at rank 6), `BW_n1_regions` (per-region axes from
the same 3-column embedding, `rec_states=regions`). Scorer `rna_state_score.py <rundir> 1.464 G10cov_B3_{1,2,3}.mrc`
(per-state axis positions in sd, deviation energy/localisation/even-odd, per-state scale on the mean). Note 7.17;
volumes in `~/Downloads/rnasep_regions_20260920/blockfit_{1,2}pc_combined/`.

**Why:** the user does not want per-axis slices ("fits that combine all the regions", "one combined embedding to
cluster from") and judged the rank-30 varimax axes unconvincing. Toy 2 showed the GMM placement needs one column per
motion and few noise columns; 1-2 PCs per region gives exactly that on the RNA.

**How to apply:** for whole-molecule RNase P states use BW_n1/BW_n2-style fits (or restrict the placement to the
leading block columns of a higher-rank fit, not yet implemented). Tail clusters (state 1 of BW_n1 at -2.4 sd, 42k
particles; state 8 of BW_n2) have lower greyscale (0.76-0.88) and are genuine low-amplitude tails, not artefacts.
Next: pose refinement against personalised references with the weights frozen, and more than 8 states.
Related: [[flex-pca-state-placement-needs-motion-frame]], [[toy2-hinge-swing-test-system]],
[[flex-pca-varimax-rotation-recovers-region-axes]].

**BW_n1_regions (same 3-column embedding, on-axis thirds):** 0.72 / 0.74 / 0.77 in-region, e/o 0.98, rel diff 0.10-0.12: as sharp as the 30-column labelled deliveries; the low-rank block fit is the RNase P default now (axes + clustered states from one embedding in ~45 min).
