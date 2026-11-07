---
name: flex-analysis-diffmap-never-worked
description: "the flex_analysis diffusion-map pipeline never produced useful results — don't spend effort preserving its behaviour"
metadata: 
  node_type: memory
  type: project
  originSessionId: 34cf3bd2-fff1-4dd4-bda2-2ec6d52c01a4
  modified: 2026-08-03T17:39:38.199Z
---

Stated by the user 2026-08-03: the `flex_analysis` diffusion-map program "never worked".
`flex_pca` (the projection-aware low-rank covariance path) is the live line of work.

**Why this matters:** you cannot tell this from the code or git history. `flex_analysis` looks
load-bearing — it is consumed by `nano3D` and `simple_trajectory_chunker`, it owns roughly half of
`src/main/flex`, and it has its own tests. That made it look like something to preserve carefully
during the 2026-08-03 cleanup, and behaviour-preservation arguments were made on its behalf that
the user did not actually care about.

**How to apply:** when weighing a change that touches both pipelines, optimise for `flex_pca` and
do not pay a design cost to keep `flex_analysis` bit-identical. Concrete instance: the RELION-style
shellwise rho floor (`reconstructor%floor_rho_shellwise`) was gated behind an optional `floor_rho`
argument specifically so the diffusion-map path kept its prior behaviour — the user's response was
that it did not matter. Still ask before deleting any of it: `nano3D` and the trajectory chunker
link against it. See [[flex-covariance-basis-defect]].
