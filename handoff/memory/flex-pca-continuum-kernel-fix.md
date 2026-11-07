---
name: flex-pca-continuum-kernel-fix
description: "10028 v6 state maps: the four continuum states were four copies of one map because GMM AUTO refits a tied-covariance mixture INSIDE a gap-free macro-cluster (responsibilities near-uniform: neff 39-67k of 73k, support 0.70 of all particles); fix = keep the equal-mass targets and give each an arc-length Epanechnikov kernel on the polyline through them (continuum_kernel_weights in simple_flex_pca_gmm), no knob; validation pending"
metadata:
  node_type: memory
  type: project
  originSessionId: 7e382bce-970c-4cdb-bc32-4a1a1c1801a9
  modified: 2026-09-08T18:45:04.758Z
---

In the default state path (state_axis=0, GMM on, GMM_AUTO on) the delivered weights are the
mixture RESPONSIBILITIES from `gmm_auto_state_weights`, not the Epanechnikov kernel whose
diagnostics the log prints first. A macro-cluster with several seats is a continuum by
construction (its components merged because `two_gauss_unimodal` found no density gap), so a
mixture refit inside it fits overlapping Gaussians: on 10028 v6, states 1-4 (15-21k particles each)
had neff 39k-67k and 70% of ALL particles in support = the whole 73k cluster, hence identical maps;
the ratchet showed only in the isolated 25k state (median distance 41) and the 6.6k island (63).
The user (2026-09-08): "the large kernel width is pretty concerning" -> go ahead with the fix.

**Fix (box branch flex-pca-resident 0387b433c):** in `gmm_auto_state_weights`, `bm == 1` macros
are unchanged; `bm >= 2` macros call `continuum_kernel_weights`: project the cluster's particles
onto the polyline through its equal-mass targets in the standardised metric
(`project_onto_target_polyline`), kernel per target on the arc-length coordinate with support =
half the target spacing (`ispace = 2*nm/bm`, the ordered-curve rule that fixed the RNA run),
`kernel_weights_at_bandwidth` with min_neff growth, labels = argmax (0 outside every support; the
outer routine now guards `lab_m >= 1`). The two helpers moved from simple_flex_pca_weights to
simple_flex_pca_util (weights re-exports them). `gmm_auto_state_weights` gained a `min_neff` arg.

**Expected on 10028:** continuum states' support ~half the cluster, neff ~20-25k, four distinct
sections along the cluster axis; states 5-6 untouched. If the user wants narrower still, the one
constant is `ispace` (nm/bm = one share). Test = states-only rerun from the v3 embedding
(`infile=flex_pca_embedding.bin`, `~/queue10.sh`, output paired_rankcut_br360_v7_kernel).
Related: [[flex-pca-gmm-undoes-equal-mass-placement]], [[flex-pca-states-stratify-by-neff]].
