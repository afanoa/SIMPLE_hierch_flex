---
name: pose-tangent-deflation-verdict-2026-11-07
description: "VERDICT lever 5: deflating the mean's rigid tangents from the flex basis fails; only per-particle E-step profiling removes pose cleanly; ported as fork wip 5f80720e0, but it costs qtoyflex closed over 5 seeds (.315 vs .441), so DEFAULT OFF since 4afe93b80 (user 2026-11-08); SIMPLE_FLEX_POSE_PROFILE=yes turns it on"
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

- **PORTED 2026-11-07 afternoon (user: "Port it").** The profile is in `embed_latents_with_contrast`, default on;
  `SIMPLE_FLEX_POSE_PROFILE=no` turns it off. Patch `handoff/pose-tangent/pose_profile_f3eab71ce.patch` (on hetero-tree-wip
  f3eab71ce); the flex-precision patch must be STACKED on it (plain git apply conflicts on adjacent imports). Not committed: the
  coordinator commits it after the gates.
  - Gates, same binary, on vs off (build_poseprof, fast gate 15/15): nulls 1 leaf (identical partitions); r131 corner kept on
    both seeds (on: s01 1.01/.52 and 1.00/.54); qtoy2 NMI .780 -> .811; cost unchanged (99 vs 96 core-h).
  - flex_pca's own IgM 5 % coordinates: shift-error R^2 .003 on vs .35 off. Off reproduces the old fflex at R^2 .998-1.000.

- **qtoyflex A/B (2026-11-08):** closed F1 on vs off, seed 1 .159 vs .474, seed 2 .413 vs .433. The site-chain session's .159 was
  this profile at seed 1.
  - Seed 1: the profile perturbs the root's flex states (this latent has no pose error), so the pole rule keeps one site. The root
    splits K=2 instead of 3, and r1's closed split falls to the effect floor (P/N 1.45 < 2).
  - The closed signal itself is intact (0.6 % in the tangent span; 4.4 % per particle removed; LDA .746 -> .745).
  - 5 seeds (2026-11-08): closed F1 on .315 +- .063 vs off .441 +- .012; holo .586 vs .601, lower in every seed; K=3 root on 0/5
    vs off 2/5. The profile thins the root's site layer (fewer kept poles), so the cost is systematic.
  - A 2-SE pose-share gate fails (qtoyflex R^2 .003 is ~14 SE); the truth-free evidence share separates cleanly (IgM .75-.87,
    elsewhere <= 0).
  - DECIDED 2026-11-08 (user: "Yes, default off"): fork wip `4afe93b80` makes the profile off by default;
    `SIMPLE_FLEX_POSE_PROFILE=yes` turns it on. Gate arms compared with earlier pp_on runs must set the key.

**Why:** settles lever 5. The port is hygiene (no IgM arm gain); it passed the gates but costs qtoyflex closed at one seed.
**How to apply:**
- Never retry basis-level tangent deflation (flex2D-style or weighted).
- Compare like with like: A/B the profile with the env toggle in one binary.
- For the core-patch shift splits, profile the site reference's in-plane shift derivatives inside the site classifier. Tested 2026-11-07: a flat profile works, a 1 px prior does not ([[site-shift-profile-verdict-2026-11-07]]).

Related: [[igm-fixture-xval-stageb-gapw-2026-11-06]], [[igm-focused-flex-mask-is-bottleneck-2026-11-07]],
[[e1-rigid-tangent-subspace]], [[site-mask-rework-direction-2026-11-07]].
