---
name: flex-pca-placement-is-deterministic
description: "flex_pca state placement has NO randomness; run-to-run variance enters through z via schedule(dynamic) accumulators, so the handoff's placement-restart plan targets the wrong stage"
metadata: 
  node_type: memory
  type: project
  originSessionId: 0e59e2de-7b03-46a6-b41f-0c4b783e5d34
  modified: 2026-08-06T19:42:16.568Z
---

Audited 2026-08-06, statically, across the whole placement path. **No RNG is reachable from any
state-placement routine.** `diffusion_kcenter_targets` picks graph nodes by fixed stride
(`simple_flex_pca_model.f90:2642`), seeds its orthogonal iteration from a deterministic sine
(`:2710`), and runs plain greedy k-center (`:2756`). `kmeans_latent_targets` seeds from the particle
nearest the latent mean then farthest-point (`:3108`). The GMM initialises at those targets and
respawns on a deterministic mean-distance test (`:1659`). Identical `z` gives identical states.

Variance enters upstream, in `z`: `!$omp parallel do schedule(dynamic)` with per-thread accumulators
— `Sbb_thr`/`SG_thr` (`simple_flex_pca_columns.f90:2764, :2868`), `hfpw_thr`/`cvpw_thr`/`gdi_thr`
(`:2779, :2783, :2862`), probe `gam_thr` (`:3887, :3929`). Fold order is fixed; the PARTITION of
particles across threads is not. Downstream choices are strict-inequality argmaxes.

**Why it matters:** `flex_handoff_2026-08-05.md` open problem 1 records ARI swinging 0.79–0.96 "at
equal embedding quality" and proposes multiple placement restarts scored by ICL. That diagnosis is
wrong — placement has no variance to remove. The variance is in `z`.

**How to apply:** to reduce run-to-run variance, change the accumulation (static `schedule`, or an
ordered reduction), not the placement. And because flex_pca is NOT bit-reproducible, a
same-configuration repeat run is not a valid null — run any comparison twice before believing a
difference. Contradicts the older assumption in [[flex-pca-speedup-verified]] that the pipeline is
deterministic enough for a repeat-run control.

Full write-up: `doc/for_developers/ideas/flex_pca_swing_regression_audit_2026-08-06.md`.
