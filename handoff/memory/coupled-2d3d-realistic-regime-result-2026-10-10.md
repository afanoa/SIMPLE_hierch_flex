---
name: coupled-2d3d-realistic-regime-result-2026-10-10
description: "toySE at SNR 0.03/0.01 with refined (biased) poses: 3D flex latent decays to a blob and loses the 10% state (recall 0.01-0.04); the pose-mean 2D class model still splits on state with refined poses (NMI 0.46/0.34, view 0.02) and pooling on its subclasses recovers the minority (recall 0.55/0.53); at SNR 0.01 one class's 2D latent beats the 3D latent per particle (0.59 vs 0.56)"
metadata:
  node_type: memory
  type: project
  originSessionId: a6c2687a-b3cd-4070-9eaa-6163aa554f8d
  modified: 2026-10-10T05:13:50.079Z
---

2026-10-10 01:00, local Mac (no VPN). Work dirs ~/toyse_local (SNR 0.1, true poses), ~/toyse_snr03, ~/toyse_snr01
(refined poses: 3 prob refine3D its at lp 8 from the truth, objfun euclid with the sigma2 file pre-placed in cwd and
1_refine3D/, mkdir=no). Chain `run_toyse_local_v2.sh`; readouts `score_toy.sh <W>`; note
doc/implementation_notes/coupled_2d3d_flex_toy_findings_2026_10_09.md; figures doc/implementation_notes/coupled_2d3d_figs/.

| | SNR 0.1 true | SNR 0.03 refined | SNR 0.01 refined |
|---|---|---|---|
| 3D latent 15-NN purity (chance 0.30) | 0.99 | 0.88 | 0.64 |
| particle k-means state-13 recall | 0.87 | 0.01 | 0.04 |
| pose-mean class model (refined poses) state/view NMI | - | 0.46/0.02 | 0.34/0.01 |
| pooled on pose-mean subclasses, state-13 recall | - | 0.55 | 0.53 |
| within-class per-particle purity 2D vs 3D | 0.73 vs 0.95 | 0.69 vs 0.76 | 0.59 vs 0.56 |

Refined poses carry the coherent per-state bias (state 13 at 7 deg vs 3 deg for the majority). The user's RNase P
observation (2D richer than 3D) reproduces only at realistic SNR + biased poses; at SNR 0.1 / true poses 3D wins.
**How to apply:** test heterogeneity claims at SNR <= 0.03 with refined poses, never at true poses; the pose mean is
the lever that makes 2D subclasses conformational; the minority state comes back through subclass pooling, not through
the per-particle 3D latent. Traps: editing a running bash chain corrupts it (write a v2 file); numpy prototypes
without OMP_NUM_THREADS pushed the load to 70; flex_pca writes state labels into ptcl3D even with rec_states=no
(selection states=1,2,3,4 merges them back); sub_prec arm of pool_place is defective. Pending: 27-state nuisance toy
~/toy27_hard (STATES=all27 NUISANCE=1 PERTURB=10 NCLS2D=100), the hardest local test.
Related: [[coupled-2d3d-pose-mean-finding-2026-10-09]], [[flex-placement-not-latent-loses-small-states]].

**Per-particle embedding (02:00):** subclass-TIED mixture on the 3D latent (isotropic noise s2/n_g on subclass means;
`tied_mixture.py`) gives posterior-mean embeddings with 15-NN purity 0.72 (SNR 0.01) / 0.91 (0.03) vs raw 0.63 / 0.87
and accuracy 0.60 / 0.67; the image E-step from the pooled seed with consensus + Wiener-shrunk deviation references
(`ptcl_posterior.py`, cross-half, fixed fractions) climbs 0.52 -> 0.61 -> 0.68 at SNR 0.01 with state-13 recall 0.90;
plain per-state references drain a state in one step (same as the 2026-10-05 soft-EM toy). The cached per-particle
precision matrices are NOT usable as measurement covariances (anisotropic tied mixture collapses; pool_place sub_prec
defective). UMAPs: coupled_2d3d_figs/umap_toyse_snr{001,03}_refined_posterior.png.

**Hard toy (03:00-04:00, ~/toy27_hard: 27 states, nuisance, 10% junk, 10-deg perturbed start, SNR 0.02):** the 3D
latent orders all three hinge angles (ridge R2 0.7-0.8, Spearman 0.9) and isolates junk (0.92); the per-class 2D latent
carries nothing (within-class angle R2 0.20 vs 0.77; purity at chance) even with mask + background/contrast columns:
the per-class basis from ~400 particles is too noisy, the shared 3D basis is not. Pooling subclasses destroys continuum
ordering (pooled/tied readouts R2 0.02-0.35). **Reading:** 3D latent = right per-particle object when the basis is good;
2D subclass pooling = the fix for lost minority states under biased poses; a 3-axis continuum is a disc in UMAP by
nature. Coupled model = shared 3D basis + class nuisance columns + subclass-tied prior (+ pose mean in the 2D step).
Next real-data test: pooled/tied placement on an existing RNase P flex embedding with the existing cls_expansion labels.
