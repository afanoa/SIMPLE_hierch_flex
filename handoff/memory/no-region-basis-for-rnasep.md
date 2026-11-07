---
name: no-region-basis-for-rnasep
description: STANDING (user, 2026-10-14): no hand-defined region / multibody (S4 four-arm label) basis replacing the whole-molecule RNase P latent; BUT data-derived FOCUSED flex (pcg_mskfile = a node or site variation mask) is wanted (user 2026-10-21)
metadata:
  type: feedback
---

On 2026-10-14 I tested the S4 four-arm region basis in the basis-span diagnostic because the 2026-10-08 memory said
"use multibody (region labels) for RNase P flex fits". The user: "i dont want this region basis. how many sessions
must I indicate that in". That memory's advice is withdrawn.

**Why:** the user's direction for RNase P is the whole-molecule basis plus the 2D-3D coupling; the region split is a
partition of a continuum (covariance ratio 1.3-1.5) and the user has ruled it out repeatedly. (The test also
showed it explains nothing of the sub-class deviations: 0.025 vs 0.017 null, automask 0.011 vs 0.010; the sphere
basis explains 0.55 vs 0.30.)

**How to apply:** for RNase P 3D modelling use the sphere-support flex_pca basis (`full692k/1_flex_pca`, lp 10
variant `full692k_lp10`) and the map-level models from the 2D sub-classes (tracked states, sliding-window
trajectory). Do not bring up pcg_mskfile / region blocks / multibody for this molecule unless the user asks. See
[[rnasep-s4-four-arm-support-test-2026-10-08]] (superseded on this point), [[rnasep-sync-od-coordinate-2026-10-14]].

2026-10-21 nuance: the user asked for "the arm-focused flex at every node" (rflex9+: flex_pca pcg_mskfile = the node's
data-derived variation mask) and, for the factorial site layer, "shouldn't we have a focused flex". So a focused flex
with a DATA-DERIVED mask (one basis per node or per site) is approved; what stays rejected is the hand-defined
multi-region basis that replaces the whole-molecule latent. Don't over-apply this note.
