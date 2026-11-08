---
name: rnasep-noise-spectrum-model-verdict-2026-11-08
description: "VERDICT 2026-11-08: RNase P's per-particle misfit is each particle's NOISE SPECTRUM (ice, denoise_project's whole-band normalisation); a 6-band empirical-Bayes noise variance (nvarB) removes the solvent floor and lifts sparse PPCA's closed AUROC (r131 span .69 -> .75-.77, 692k-50k .79 -> .85); knob-free, collapses on IgM; port into flex's whitening started"
metadata:
  node_type: memory
  type: project
  originSessionId: cc23038c-fea9-4b5f-b321-d58eb2942350
  modified: 2026-11-08T07:21:09.619Z
---

rnasep-misfit session; batch2 `hetero_tree_20261028/handoff/rnasep-misfit/{VERDICT,STATUS}.md`, runs `math_reworks/rnasep_misfit/`.
- The misfit: residual power over 117-15 A has a true SD of 29 % (reliability .96) on r131 and 692k. The mean, nuisances and sites
  remove ~2 % of it, and it carries no closed signal. It is low-frequency, CTF-shaped and clustered by FoilHole / micrograph (ice),
  anti-correlated with 6-3 A power (denoise_project normalises each image by its whole-band sd). The 2D registration is NOT it.
- The model (nvarB): an inverse-gamma prior per band (6 bands) fitted by marginal likelihood; each particle's posterior-mean variance
  multiplies sigma2 in the whitening and in the fission draw. Knob-free; on IgM the prior CV is <= .013, so it does nothing there.
- Scores, base -> nvarB:
  - r131 floor: balls outside the molecule 76 -> 0.4 nats (z 0.3); rim controls 245 -> 29. A scalar per particle fixes 20 %; a band
    scale and a registration correction nothing.
  - Latent: log-likelihood vs noise r .61 -> .12; closed at the fixed native basis unchanged (.706 -> .710).
  - Sparse pursuit WITH a per-round control check: r131 best block .640 -> .70-.74, span .694 -> .75-.77 (the same site every run);
    692k 50k subsample .737 -> .791 / .794 -> .853, beyond the native run on all 692k (.802 / .849). IgM arms within .003.
  - The per-round check is REQUIRED: without the floor, the as-written stop runs away.
  - Shared tree at the 692k subsample: the noise share of the root split goes .20 -> .07. The root then splits along amplitude
    (.38): a 10 % low-amplitude, pose-unstable side. Closed is not at level 1 (AUROC .495).

**Why:** this is the observation-model fix the shared tree model needed on real data ([[real-data-fission-needs-misfit-model-2026-11-08]]),
and it helps sparse PPCA on its own.
**How to apply:**
- Real-data flex and sparse fits use a per-particle band-wise sigma2. Port: flex already whitens per particle
  (`build%esig%sigma2_noise(shell, iptcl)`), so only the values change: one residual pass at the node mean + 6 two-parameter fits.
  The stkind sigma2 groups on the RNase P stacks are the 10 denoise parts, not micrographs.
- Do not call the 692k low-amplitude side junk: the user judged that group real RNA conformations ([[rnasep-junk-screen-off-2026-10-23]]).
  Read it by honest maps and by growing below it.
Related: [[tree-latent-shared-basis-verdict-2026-11-07]], [[shared-tree-model-rnasep-direction-2026-11-07]], [[band-must-be-data-derived]].
