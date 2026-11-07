---
name: rnasep-bins80-within-view-motion-2026-10-15
description: "first positive direction-bin result on RNase P (corrected bins80, 2026-10-15 05:00): within one projection bin the cls_expansion sub-classes differ by a reproducible, arm-localised 20-40 A difference that survives rigid in-plane registration and pose-tangent deflation (even/odd 0.76-0.94); the pooled 3D latent grid at consensus poses stays diffuse"
metadata:
  node_type: memory
  type: project
  originSessionId: 5b04614c-23e0-4650-8600-03b464a7f576
  modified: 2026-10-15T08:51:50.225Z
---

Corrected bins80 run (row-ordered bin-frame stacks, job 53878700, 11 min): 100 sharp distinct bin averages, 1,599
conditional sub-classes. Reconstruction-free readout `bins2d_dipole.py` (scripts/cls_expansion_proto/subclass_pool):
per bin, upper minus lower quartile of sub-classes along the within-bin latent z1, band 100-15 A, even/odd from the
even/odd sub-class stacks. Median of the 20 biggest bins: raw rms 0.66 of the average, even/odd 0.89, lobes 22 A
apart; after rigid registration (+-15 deg / 5 px) and deflation of the in-plane pose tangents still rms 0.39,
even/odd 0.76, lobes 29 A, share 0.28; big bins 0.8-0.94 with red/blue pairs on the arms (bin 22 both arms 30 A,
bin 95 upper arm 43 A, bin 69 36 A). A third of the raw difference is rigid in-plane disagreement (the within-bin
latent's first mode is partly residual in-plane rotation).

The pooled lp10 flex latent per bin sub-class is reproducible (within-bin PC1/PC2 even/odd 0.976/0.960) but its 3 x 3
grid reconstructed at the consensus poses is the usual diffuse speckle (lobes 1-6 %), so 3D at consensus poses is
not the readout; the synchronised 2D coordinate (sync windows bsyn1-3, grid bsg, 2D ladders sync_ladder.py) is.
Figures: doc/implementation_notes/rnasep_coupled_figs/viz_rnasep_bins80_{cavgs,subclasses_early,2d_dipoles*,pooled_grid}.png;
notes 21.4-21.6.

**Update 05:25:** the cross-view tests (near_pair_consistency.py, near_pair_subspace.py) are null: neighbouring views
(< 20 deg) do not share the leading within-bin difference (0.03 signed, 0.09 free sign) nor any deviation mode
(top canonical correlation 0.615 near vs 0.605 far); the bsyn windows and bsg grid at refined poses give even/odd
0.05-0.24 with 0-2 % lobes. Not defocus (corr -0.01), not tilt (E_view 0.09). The within-bin even/odd of 0.76-0.9 is
consistent with selection bias (quartiles selected on the images themselves). Variants with fewer sub-classes and a
coarser band (4 @ 30 A, 8 @ 25 A, 4 @ 40 A) were launched to see whether a real large motion can dominate.

**Calibration (05:45):** the same tests on the toy bins (toyL3/BINS/toy, real hinge motions, 400 particles per bin)
give near-pair leading-mode |corr| 0.50 and subspace overlap 0.78 near vs 0.53 far, so the tests discriminate; RNase P
stays at 0.09-0.12 for 16 @ 20 A, 4 @ 30 A, 8 @ 25 A, 4 @ 40 A and for 50 bins (8 @ 25 A): null at every bin size and band. In-plane misalignment is 3 % of members
(subclass_misalignment.py). Scripts: near_pair_consistency.py (--key pca needs no latent dumps),
near_pair_subspace.py, subclass_misalignment.py, run_bins80_variant.sbatch (DIR/NSUB/LP2D/ORITAB).

**Why:** the only immune metrics are cross-view ones; the within-bin reproducibility fooled me for an hour.

**How to apply:** never report even/odd reproducibility of a difference between image-selected subsets as evidence;
use cross-view consistency (near-pair tests) or an independent-noise band with enough SNR. Judge direction-bin
variants by the near-pair leading-mode |corr| (0.09 baseline). See [[rnasep-692k-project-rows-are-a-permutation-of-stack-order]],
[[rnasep-bins-overnight-2026-10-15]], [[toyl-direction-bin-scheme-results-2026-10-14]].
