---
name: rnasep-final-poses-make-states-separable
description: "2026-09-28 RNase P: refine3D_states FINAL labels go Mahal^2 0.19 -> 3.74 when the latent is re-embedded on the final poses; CIRCULARITY CONFIRMED on a homogeneous toy 2026-10-04 (null labels 0.71 -> 3.15), so this is not evidence of real states"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T06:23:52.060Z
---

Test (job 53381910, `/mnt/beegfs/elmlund/afan/rnasep_soft/FF`, script `~/rnasep_flexfinal.sbatch`): flex_pca with the
refine3D_states flex-init settings (npreimages=5 gridding min_state_frac=0.1 sigma_est=global mskdiam=200, rec_states=no)
on `rnasep_states/25_refine3D_states` final poses (labels merged to 1, consensus from final poses as vol1), scored by an
LDA on the latent (state means, pooled covariance) against the refine3D_states FINAL labels.

| latent | labels | pair Mahal^2 min/med/max | argmax acc (chance) |
|---|---|---|---|
| seed (consensus poses) | final | 0.02/0.19/0.29 | 0.36 (0.23) |
| seed | own flex | 0.38/0.50/1.17 | 0.41 (0.21) |
| final poses | final | 1.24/3.74/5.27 | 0.59 (0.23) |
| final poses | own flex | 0.06/0.16/0.29 | 0.34 (0.20) |

**Reading:** pose refinement per state makes the domain states visible to the linear latent (~20x Mahal^2), the user's
"couple refinement to heterogeneity" hypothesis. Two open caveats: (1) circularity — final poses were refined against
state-specific maps, so labels imprint on poses (reference bias) and the separation may be partly induced; needs a
reproducibility control (independent stochastic-init refine3D_states run) before believing it; (2) flex's own placement
(diffusion k-center + population floor) does NOT recover the separable states even in the new latent, so placement, not
only the latent, is a weak link.

Also measured: soft state posteriors from the seed partition are nearly flat on RNase P (median max posterior 0.33,
95-99.6% of particles >= 0.05 in every state), so the per-state soft refinements S1-S5 (`rnasep_soft/S*`, refine3D_auto
flex_weights=yes flex_state=s, 25k/iter, maxits 50) are close to reweighted consensus refinements; the user chose to run
them anyway to see whether per-state poses alone drive domain separation. Related: [[rnasep-flex-seeds-vs-refine3d-finals]],
[[rnasep-flex-per-state-plan]].

**Circularity test, 2026-10-04 ([[pose-state-coupling-toy-verdict]]):**
- Setup: the same refine3D_states invocation (`nstates=4 flex=yes mskdiam=200`, `~/SIMPLE_flexref6`), run on a
  HOMOGENEOUS toy (toy3 at theta = 0 only, 40k, SNR 0.02). Batch2 `~/toyPS/r3s.sbatch`, dirs
  `toyPS/t3null_r3s` and `toyPS/t3_r3s`.
- On the null it returned four "states" (7.5k / 4.6k / 8.4k / 19.4k).
- Their separability rose from Mahal2 median 0.71 in the latent at the starting poses to 3.15 in the latent
  re-embedded at its final poses.
- On toy3 with a real swing the rise is the same size: 5.5 -> 7.5.
- **Reading:** per-state pose refinement imprints any partition on the poses. The RNase P 0.19 -> 3.74 is what
  imprinting alone produces, and caveat (1) above is the explanation.
- On toy3 refine3D_states also assigned states WORSE than plain latent bins: theta variance explained by 4
  states 0.774 vs 0.841 (oracle/true poses 0.882).
