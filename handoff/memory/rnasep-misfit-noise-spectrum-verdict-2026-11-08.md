---
name: rnasep-misfit-noise-spectrum-verdict-2026-11-08
description: "VERDICT: RNase P per-particle misfit = each particle's in-band NOISE SPECTRUM (sd 29 %, reliability .96, FoilHole/ice, CTF-shaped), not registration; per-particle band-wise EB noise variance kills the solvent floor and lets sparse PPCA read closed (r131 block .70-.74, span .75-.77 vs .66/.71); IgM unchanged"
metadata:
  node_type: memory
  type: project
  originSessionId: 83276503-c95a-4f22-abab-30791a97562c
  modified: 2026-11-08T07:56:49.764Z
---

rnasep-misfit session, 2026-11-07/08 (batch2 `hetero_tree_20261028/handoff/rnasep-misfit/VERDICT.md`, STATUS.md; work dir
`math_reworks/rnasep_misfit/`).

- What the misfit is: against one white noise variance per shell, an r131 / 692k particle's residual power in 117-15 A has true sd
  29 % (split-half reliability .96), untouched by the mean + nuisances + sites, not closed (AUROC .48-.51). Low-frequency structural
  noise: tracks the image's in-band power, anti-correlated with its 6-3 A power (denoise_project normalised each image by its
  whole-band background sd), FoilHole / micrograph clustered (ice; micrograph particle count, not neighbours in the box), CTF-shaped.
  The 2D registration (rotation, shift, sub-pixel) is NOT involved (|r| <= .08; phase-ramp shift, KB rotation, no padded corners).
- It is what made "any local mode gains ~100 nats" and tree-latent's 692k misfit split (child 1: +24 % low-frequency noise, child 2:
  amplitude .20 = signal-free particles).
- Fix (`nvarB`, `mf_obs.py`): per-particle noise variance per band (6 bands) = posterior mean of an inverse-gamma prior fitted per band
  by marginal likelihood to the whitened residual sums. Knob-free; collapses on stationary noise (IgM prior CV <= .013).
  - Floor: balls outside the molecule 76 -> 0.4 nats (z 0.3); rim controls 245 -> 29 (nvarBC, + a CTF^2 noise slope: 16).
  - Latent no longer tracks the noise (ll r .61 -> .12).
  - With the per-round grown-block check (REQUIRED once the floor is gone: else the pursuit runs away), r131's pursuit stops at
    7-9 blocks and reads closed: best block .743 / .696+.649 / .730 over 3 runs, span .761 / .766 / .749, tail share .17-.21.
    Old model: native best .659, span .711; Python base .640 / .694. The closed block sits at ~(40, -33, 7) A, 50-60 A out.
  - IgM 10 % arms within .003. 692k 50k-subsample pursuit: best block .737 -> .791, span .794 -> .853 (native needed all 692k).
  - Shared tree (tree-latent's rule): at r131 one leaf on every basis / model (honest map .81 / .62; K=2 XD gain .127 -> .0095:
    misfit gone, closed is d' ~.9, too weak for a mixture). At the 692k subsample the root still splits, now off a 10 %
    pose-unstable, low-amplitude population (not noise). Do NOT call it junk: the user judged the low-amplitude group at the
    692k root real RNA conformations ([[rnasep-junk-screen-off-2026-10-23]]), so read it by honest maps and growth below it.
    Next tree lever: separability-ranked proposals.
- Port: flex already whitens per particle (`esig%sigma2_noise(shell, iptcl)`, filled in `load_and_validate_sigma`); only the values
  change (one residual pass + 6 two-parameter fits). TRAP: on these stacks the sigma2 groups (stkind) are the 10 denoise parts.

**Why:** the RNase P blocker for the shared tree model and sparse sites was per-particle misfit; this is its cause and fix.
**How to apply:** any real-data latent / site / fission work on the RNase P stacks needs per-particle band-wise noise first; draw
fission noise at sig2 x v_ib. Related: [[shared-tree-model-rnasep-direction-2026-11-07]], [[band-must-be-data-derived]], [[tree-latent-verdict-2026-11-07]], [[extremes-fission-verdict-2026-11-07]],
[[rnasep-denoised-stacks-are-phase-flipped]].
