---
name: sigma2-bootstrap-group-empty-half
description: "refine3D on a project without a registered sigma2 store dies in the capped bootstrap with 'canonical sigma2 group has an empty even/odd half' under the default per-stack grouping; pass sigma_est=global (as the flex recipes do)"
metadata:
  node_type: memory
  type: project
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-10-25T21:25:05.762Z
---

2026-09-21, master 08bd48202+: `simple_exec prg=refine3D` on the 10180 4_abinitio3D project copy (no
`sigma2_state` registered) ran `SIGMA2 INIT: rebuilding canonical state`, sampled 12500/12500 even/odd of
23375/23408 (CAP 25000), then aborted: `canonical sigma2 group has an empty even/odd half;
simple_commanders_euclid_distr.f90 line 223`. Cause: with the default per-stack grouping (`sigma_est=group`)
the capped bootstrap subsample leaves some stack with no particle in one half, and
`sigma2_state_reduce_groups` treats an empty half as an error rather than falling back.

UPDATE 2026-10-25: `sigma_est` is no longer a command-line key in release-4 builds (SIMPLE_proc aea5b1846 rejects it as
"unknown command-line argument"); global sigma estimation is the default since 9e0f6e237 (2026-10-01). Do not pass it.
**How to apply (pre-release-4 builds only):** add `sigma_est=global` to refine3D on such projects (the flex_pca recipes already do), or
seed the store first. Worth a fix in the bootstrap sampler (per-group coverage or a global fallback);
not fixed as of this note. Related: [[flex-weights-store-2026-09-21]].
