---
name: flex-pca-must-be-knob-free
description: "Standing constraint — flex_pca must work agnostically on RNA continuum AND synthetic discrete data with no data-specific knobs; state_axis/min_neff/nbins must become internal and auto, aiming at real-time heterogeneity"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 126275fa-8460-444b-8ba7-1e0140968435
  modified: 2026-08-13T16:53:06.093Z
---

Stated 2026-08-13, while diagnosing the RNA continuum failure (run 29/30): the user wants ONE
method that works on the afm RNA data and all synthetic sets (IgG, Ribosembly, Tomotwin) without
per-dataset tuning. `state_axis`, `min_neff`, `nbins` (and by extension bandwidths, K, ceiling)
"shouldn't be user facing and should be automatically determined based on the data."

**Why:** the target is real-time heterogeneity analysis (stream/refinement-integrated). A method
tuned per dataset by a human cannot run during collection. Corollary the user did not state but
that follows: auto-tuning must be CHEAP — it has to live in latent space (seconds), not in volume
space. The current nbins CV selects bandwidth by ~10 distributed trial reconstruction rounds
(~15-30 min at 932k) — right idea, wrong cost tier for the end state; it needs a latent-space
surrogate (held-out likelihood / split-half stability).

**How to apply:** treat hand-set placement knobs in runs as DIAGNOSIS ONLY (run 30 style), never
as the deliverable. Design decisions should route toward: equal-mass (quantile) target placement
(unifies discrete + continuum), cross-halfset component reliability for auto axis/rank selection,
neff-targeted bandwidths with internal adaptation on frac_in_support, over-provision + complete-
linkage merge as the K estimator, shrinkage always on. `preimage_auto` is a stepping stone; the
end state is that auto IS the default and the flag disappears. See
[[afm-rna-flex-pca-needs-continuum-readout]], [[flex-pca-signal-ranking-is-blocked]] (the split-half
trap that any auto-ranking must avoid).
