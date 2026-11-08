---
name: rnasep-misfit-noise-spectrum-verdict-2026-11-08
description: "VERDICT: RNase P per-particle misfit = each particle's in-band NOISE SPECTRUM (sd 29 %, reliability .96, FoilHole/ice, CTF-shaped), not registration; per-particle band-wise EB noise variance kills the solvent floor and lets sparse PPCA read closed (r131 block .70-.74, span .75-.77 vs .66/.71); IgM unchanged; committed default-off as 03506060d; r131 corner 3/4 seeds on vs pp_off 4/4, closed purity .099 vs .079"
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
- Native port (2026-11-08, NOT committed; patch `handoff/rnasep-misfit/noisevar_5b17de846.patch`, build
  `wip_builds/noisevar_5b17de8`, worktree `/Users/atifao/SIMPLE_noisevar`): `esig%sigma2_noise` is (shell, eo HALF), not per
  particle, so the port scales each particle's observation planes per band in `planes_batch_load` (`flex_plane_store%whiten_bands`);
  `simple_flex_noise_var` estimates in prepare, workers adopt `flex_pca_noise_var.txt`; `SIMPLE_FLEX_NOISE_VAR=no` reproduces pp_off
  exactly. Native r131: floor 109.5 -> 9.4, span .711 -> .765, best block .659 -> .681; IgM identical; qtoyflex within .001.
  COMMITTED default-off (`SIMPLE_FLEX_NOISE_VAR=yes` turns it on) on fork `hetero-tree-wip` as 03506060d: band edges are fractions
  of kb (today's rings at kb 15), forced-team tester `flex_noise_var`, workers follow the master's file. r131 corner over 4 seeds:
  pp_off 4/4, noise variance 3/4 (seed 1 .947, seed 4 .953), but its best depth-1 group holds more closed reference rows on every
  seed (purity .099 vs .079); qtoy2 seed 2 NMI .778 (seed 1's 43/49 merge was a draw). Default-on = the user's decision. Native = Python up to flex's own per-image background normalisation (r .994-.999 after it); flex floors the
  band limit (ring 15 at lp 15), the Python framework rounds it (16). On toys the node-mean residual takes a few % of strong state
  contrast into v (qtoy2 holo AUROC of v .61). TRAP: on these stacks the sigma2 groups (stkind) are the 10 denoise parts.
- Alternatives scored: micrograph / FoilHole groups and the background-window estimate (6 bands) re-learn the closed block at
  .72-.73 too but leave 10-25x nvarB's outside floor; 2 residual bands fail. See [[latent-framework-lp-sets-lattice-overwrites-cache]].

- COMMITTED 2026-11-08 as fork wip `03506060d`, OFF by default (`SIMPLE_FLEX_NOISE_VAR=yes` turns it on); own build + fast gate
  15/15. Default-on waits for the user: the r131 corner read .947 on seed 1 (bar .95), within pp_off's own seed spread, and that
  own-data readout is largely selection imprint.
- 4-seed r131 table (2026-11-08 morning): own-data corner nv 3/4 vs pp_off 4/4 (seed 1 .947), but the independent closed-reference
  purity of the best depth-1 group is higher with nv on every seed (mean .099 vs .079; recall .60-.80 vs .54-.62); qtoy2 seed 2
  equal (.778 vs .780). Coordinator recommended default on; the user's call.
**Why:** the RNase P blocker for the shared tree model and sparse sites was per-particle misfit; this is its cause and fix.
**How to apply:** any real-data latent / site / fission work on the RNase P stacks needs per-particle band-wise noise first; draw
fission noise at sig2 x v_ib. Related: [[shared-tree-model-rnasep-direction-2026-11-07]], [[band-must-be-data-derived]], [[tree-latent-verdict-2026-11-07]], [[extremes-fission-verdict-2026-11-07]],
[[rnasep-denoised-stacks-are-phase-flipped]].
- DECIDED 2026-11-08 (user: "Yes, default on"): the noise model becomes the default on wip, with SIMPLE_FLEX_NOISE_VAR=no to turn
  it off.
  - Basis: the honest 4-seed r131 gate does not separate the arms (1/4 off vs 2/4 on). The noise model shrinks the own-data
    imprint (.04-.06 vs .11-.18) and removes the misfit floor.
  - rnasep-misfit flips it, builds (nvon_<hash>), fast-gates and smoke-tests it; then cherry-pick to wip.
- ON WIP 2026-11-08 evening as 938ed9a8d (cherry-picked from dd02eff12; build nvon_dd02eff12, fast gate 15/15, build_lf's OpenBLAS
  pair). The smoke runs with the variable unset show the master estimating in both flex steps and all 4 workers adopting.
  SIMPLE_FLEX_NOISE_VAR=no turns it off.
