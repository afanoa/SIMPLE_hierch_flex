---
name: pose-tangent-deflation-verdict-2026-11-07
description: "VERDICT lever 5: deflating the mean's rigid tangents from the flex basis: plain metric does nothing, information metric costs qtoy2 signal; only per-particle E-step profiling removes pose cleanly, but it changes nothing on the IgM failure (port held)"
metadata:
  node_type: memory
  type: project
  originSessionId: 42503210-8e5e-4af5-8adc-0b1757c44b8f
  modified: 2026-11-07T17:14:42.124Z
---

Lever 5 of guide 7.8, run 2026-11-07 (Python, `latent_states_v4` statistics, the node's own fflex basis, no refit). Write-up:
`hetero_tree_20261028/handoff/pose-tangent/VERDICT.md` on batch2 (`tables.md`, `treesplit.md` beside it).

- **The IgM roots' whole-molecule latent is pose.** Fixture start error 1.43 deg / 4 A. The leading axis is 25-30 % shift
  error; the production Fortran z's are 27-39 % (shift R^2 .30-.35). Per-particle profiling removes 75-87 % of the latent's
  evidence. The nulls' only mixture gain (K4 +.002 nats/particle) is pose.
- **qtoy2 and r131 latents are not pose-contaminated** (shift R^2 <= .01; nothing tracks stage B's corrections).
- **Basis deflation fails.**
  - Plain metric (flex2D `d1a20eff6`, and the Gram-Schmidt of flex's mean-shaped deflation set in
    `simple_flex_probe_fit_update.f90`): no effect at all.
  - Information metric (D, or its per-shell average Dr): removes the pooled axis, but costs qtoy2 closed .964 -> .944 and combo
    NMI .50 -> .46. Domain motions have rigid components.
- **Per-particle profiling (the 6 projected tangents Schur-complemented in the E-step) is clean.** Pose tracking goes to 0;
  qtoy2 is kept or slightly better; r131 closed holds.
- **Still no IgM gain.** The arms stay invisible. The tree's false IgM root splits do not follow pose error (AUROC .51-.54).
  The only pose-driven IgM split is the site classifier's core patch (AUROC .81).

**Why:** settles lever 5. The port was held: it is hygiene, not a fix for the hard test.
**How to apply:**
- Never retry basis-level tangent deflation (flex2D-style or weighted).
- If pose hygiene is wanted, port the profile into `embed_latents_with_contrast` (about 150 lines).
- For the core-patch shift splits, profile the site reference's in-plane shift derivatives inside the site classifier. Tested 2026-11-07: a flat profile works, a 1 px prior does not ([[site-shift-profile-verdict-2026-11-07]]).

Related: [[igm-fixture-xval-stageb-gapw-2026-11-06]], [[igm-focused-flex-mask-is-bottleneck-2026-11-07]],
[[e1-rigid-tangent-subspace]], [[site-mask-rework-direction-2026-11-07]].
