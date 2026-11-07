---
name: flex-pca-state-placement-needs-motion-frame
description: "2026-09-20 toy: whole-molecule GMM states (8, full latent) resolve only the dominant motion from a 12-column block fit or global fit (pairwise-difference truth 0.66) but reproduce joint theta2 x theta3 states from the rank-4 varimax frame (0.92, as distinct as the truth); the standardised metric weights every column equally so noise columns swamp weak motions -- place states on region-attributed / reliability-weighted columns, not on all columns"
metadata:
  type: project
---

Scorer: `scripts/flex_pca_verify/region_toy/toy_state_score.py <rundir> <toydir> [outdir]` -- each delivered
state vs the EXACT reference (true volumes of its particles, same weights); whole-map corr is useless (0.98
everywhere, the mean dominates); use the pairwise (state_s - state_t) vs (ref_s - ref_t) and the
deviation-from-set-mean tests. Runs: toy/T_B3_states8, T_G12_states8, T_G4vmx_states8 (npreimages=8,
pcg_mskfile=labels96_B3 / env96, SIMPLE_COV_VARIMAX=1 for the last). Note section 7.15.

**Why:** the GMM placement standardises every latent column to unit variance; 9 of 12 block-fit columns are
noise/leak columns and get the same weight as the 3 motions, so only the 12-degree motion separates
clusters. The rank-4 varimax frame has one column per motion.

**How to apply:** for whole-molecule states on RNase P in the rank-30 varimax frame, restrict the placement
to the region leaders (e.g. the 6 columns of 7.14) or weight columns by their localisation/reliability;
the on-axis per-region deliveries are unaffected. The user wants whole-molecule states, not axes
("fits that combine all the regions"). Related: [[flex-pca-varimax-rotation-recovers-region-axes]],
[[flex-pca-states-stratify-by-neff]].
