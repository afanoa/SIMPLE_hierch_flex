---
name: cls-expansion-vs-flex-pca-latent-diagnosis-2026-10-09
description: "Investigation (no code) of why cls_expansion's 2D latent separates RNase P states while flex_pca's 3D latent does not: view-conditioning + subclass pooling vs global linear-first; proposed links = subclass-pooled 3D placement, 3D covariance from per-class 2D deviations, image E-step, refine last"
metadata:
  node_type: memory
  type: project
  originSessionId: a6c2687a-b3cd-4070-9eaa-6163aa554f8d
  modified: 2026-10-10T01:20:54.235Z
---

2026-10-09, user question: cls_expansion (branch cls-expansion-solvent) gets every RNase P state in 2D, flex_pca
does not; why, and how to couple them. Investigation only.

**Diagnosis (from recorded measurements):** flex_pca is linear-first (global subspace at consensus poses, mixture
on a 99.7%-noise per-particle latent, standardised GMM, macro-cluster fold); cls_expansion is mixture-first
(abinitio2D fixes view x coarse state) with a per-class PPCA, in-plane pose tangents as nuisance, solvent prior, and
the delivered unit is a subclass of ~100-500 members = exactly the group size the IgG SNR analysis said to pool.
3D entangles view and heterogeneity operationally: view-dependent visibility smears conformations into a continuum,
pose bias is coherent per conformation, and weak-signal partitions default to view/contrast splits. The Nature paper
got states by discrete 3D classification (image-likelihood mixture EM), the same order.
**Caveats:** 2D richness is partly readout (model renders look clean; pre-prior main axis was background);
the only measured 2D->3D link (solve3D_cavgs docking) is chance-level, so 3D consistency is unproven.

**Proposed links, cheapest first:** (1) pool flex_pca z over cls_expansion subclasses and place on subclass means
(rerun the toy3 cluster-tracking test with subclasses, not classes); (2) fit the 3D covariance to per-class 2D
deviations U_v z (parent cancels, deflate cone-spread rigid tangents), Python toy first; (3) finish with the soft-EM
image E-step and refine3D_auto per state LAST; never init with refine3D_states flex=yes (sphere+gridding).
Related: [[cavg3D-states-findings-2026-10-08]], [[flex-placement-not-latent-loses-small-states]],
[[pose-state-coupling-toy-verdict]], [[flex-pca-latent-is-997-percent-noise]], [[cls-expansion-solvent-prior-2026-10-09]].
