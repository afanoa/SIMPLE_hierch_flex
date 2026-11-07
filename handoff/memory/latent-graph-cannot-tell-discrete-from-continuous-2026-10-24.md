---
name: latent-graph-cannot-tell-discrete-from-continuous-2026-10-24
description: "Phase A, the user's discrete-split vs continuous-swing direction: diffusion maps + sweep conductance on root flex/fflex MAP z show NO bottleneck even for discrete states (CAHRA open/closed, cdisc A ratios 0.92-1.06); an imposed two-peaked IgM occupancy is lost (ratio 0.97). Per-particle latent noise ~ signal; decide on averaged (binned-map) evidence: per-shell mixture residual + step share (Phase B-C); OLE's 3 series call two-state (Phase D)"
metadata:
  node_type: memory
  type: project
  originSessionId: 1ca67e82-4d15-4725-b95d-0ee5719be9f4
  modified: 2026-10-25T03:35:37.002Z
---

2026-10-24, for the coordinator (user direction: keep the tree, but at each split tell a discrete change from a continuous swing; graph
route = diffusion maps + Cheeger conductance). Tool: routes/toy55 `dmap_test.py`; results in rnasep_klem/factorial/validation/dmap_phaseA/
and factorial_site_layer_results_2026_10_21.md.
- Synthetic low-noise control works: a junk cluster or two-state switch gives a conductance ratio of 0.00 vs a Gaussian null; an
  imposed dip 0.14; a uniform swing ~0.8 (not 1, since uniform has a lower middle cut than a Gaussian).
- Real roots (k = 15 / 30): CAHRA M s1 ratio 0.92-1.04, cdisc M 0.94-1.06, IgM 0.50-0.89; every profile unimodal.
  - Diffusion coordinates track truth more weakly than raw latent components: IgM a4_th 0.56 vs fflex z2 0.72.
- An IgM fflex subsample with an imposed two-peaked a4_th occupancy (dip at 8-13 % of peak) reads ratio 0.97, one mode, KS 0.30.
- Flex's deconvolved GMM AUTO also misreads: it flags CAHRA and cdisc fflex as continuum, and IgM as 8 macro-clusters.
- PHASE B, binned maps (phaseB_maptests.py; rnasep_klem/factorial/phaseB/):
  - The step profile works qualitatively. Switches show a dominant step at the switch bin (CAHRA open 0.61, cdisc closed 0.60, cdisc
    A 0.47); IgM swings are even (0.28-0.39).
  - The end-map mixture test (bins 2-4 as a P1 + b P5) does NOT separate on real data, though it does on synthetic maps. Every series
    has a reproducible residual (z 18-65), and beta misses the true fraction (CAHRA 0.20-0.54 vs 0.02-0.24), because bins along a
    score also differ in junk share and per-particle contrast.
  - Proposed fix: amplitude-match per shell and use a moving-region focus mask before the fit.
- PHASE C (phaseC_maptests.py):
  - Per-shell amplitude matching DISTORTS true blends (synthetic clean switch 1.0 -> 2.3-3.5 x noise); don't use it before a
    mixture fit. The focus mask adds nothing; for CAHRA it is 87-89 % of the molecule (a global difference).
  - The per-shell mixture fit (a(s) P1 + b(s) P5 per Fourier shell) is blur-robust and gives a clean gap. In the root mask at
    15 A, switches are 1.7-5.0 x noise vs swings >= 11.5 (10 A: 1.3-2.7 vs >= 6.25). The focus mask stretches the switches to 9.8
    (swings >= 12.7).
  - Step share: switches 0.47-0.61 vs swings 0.28-0.39 at 15 A.
  - beta still misses the state fraction. For CAHRA it is NOT brightness (phase E: contrast does not track open, and P5 - P1 is not a
    rescaled map); for cdisc, junk varies across bins. See [[at-node-switch-swing-call-2026-10-24]].
  - The call to use: per-shell residual + step share. Check beta only on cleaned scores.
- PHASE D, OLE real data (phaseD_maptests.py; ole_12707/motion_series/switch_swing.txt; 0.1 core-h):
  - The residual puts all three series in the switch range: r1_s4 1.3-2.2, r2_s4 1.2-1.9, r_s2 2.6-3.8 (root mask, 15 A).
  - Only r1_s4's step share is switch-like (0.48); r2_s4 (0.30) and r_s2 (0.33) sit in the swing range.
  - r1_s4 and r2_s4 are one-signed occupancy changes on the C2 axis at z -58 to -65 A. They are the same feature seen from opposite
    poles (posterior Spearman -0.63).
  - r_s2 is 83 % gain at z +78 A with a minor negative lobe; its call is switch-leaning, lower confidence.
  - Blind spot: every toy swing is a large motion. A shift well below the resolution gives blend-like intermediates and reads as a
    switch, so check the sign balance of P5 - P1 (one-signed = occupancy).
**Why:** the best single latent component tracks a motion at abs(rho) 0.5-0.7, so the per-particle noise is about the signal.
Particle-level density cannot show gaps.
**How to apply:** don't build the discrete / continuous decision on particle-level latent graphs or densities. Use averaged
evidence, such as the step profile of a binned map series ([[motion-series-prototype-2026-10-24]]): one dominant step for a switch,
even steps for a swing. Related: [[discrete-splits-capture-the-span-states-that-matter]], [[igm-motion-extremes-capture-2026-10-24]].
