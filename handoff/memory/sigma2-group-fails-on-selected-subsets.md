---
name: sigma2-group-fails-on-selected-subsets
description: "flex_pca on a 2D-selected subset dies with 'canonical sigma2 group has an empty even/odd half' - per-stack sigma2 is structurally impossible there; pass sigma_est=global"
metadata:
  node_type: memory
  type: project
---

Running flex_pca (or anything euclid/OBJFUN_EUCLID) on a post-2D-selection project fails at startup with
`ERROR! canonical sigma2 group has an empty even/odd half; simple_commanders_euclid_distr.f90:223`.

**Why:** with the default `sigma_est=group`, `ngroups = max(stkind)` over ACTIVE particles
(`simple_commanders_flex_pca.f90` ensure_canonical_sigma_state, and `simple_sigma2_files.f90:60-70`),
and `sigma2_state_reduce_groups` rejects the state if `any(counts == 0)`
(`simple_sigma2_state.f90:288`). 2D selection deselects EVERY particle in some stacks, so those stack
indices survive as groups with zero members. It is not a bug: you cannot estimate a per-stack noise
spectrum for a stack with no surviving particles.

**Fix:** add `sigma_est=global` (sets `params%l_sigma_glob`, ngroups=1). Verified on EMPIAR-10180:
the 46,783-particle selected subset failed with group sigmas and ran with global
(`FLEX_PCA particles=46783`). A full-population run never hits this (all stacks populated), which is
why the 327k realign run worked with the default.

**How to apply:** any time you restrict to a selection, expect this and pass sigma_est=global. NOTE it
is a real deviation from a full-population baseline's noise model - say so when comparing results. The
alternative that preserves group sigmas is to run from the FULL project (all stacks populated) and
restrict the fit with `pindfile=`, since ngroups comes from the project's state field, not the pindfile.
See [[flex-pca-validate-inputs-before-batteries]], [[spliceosome-realign-deselected-unlocks-basis]].
