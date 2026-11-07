---
name: pooled-2d-route-sees-small-domains-synthetic-2026-10-17
description: on the LIGHT synthetic sets (12 kDa arms / 11 kDa subunits at snr 0.0012) the pooled 2D route (abinitio2D ncls 100 -> cls_expansion kernel sub-classes of 5 classes, ~100-190 particles each) carries the small-domain signal the klem boards cannot see: best arm std ratio 0.46 (klem leaves 0.74), best site occupancy 0.85/0.15 (klem boards none)
metadata:
  type: project
---

2026-10-17 19:35, route2d_arms / route2d_tethered (routes/toy55/route2d_<mode>/, synth_2droute.sbatch: fresh project
re-imported with ~/SIMPLE/build, abinitio2D ncls 100 mskdiam 200 nparts 10 -> clsexp_stage1/run_kernel_cls.sh on the 5
best classes with CLASSES exported -> the ptcl2D `class` field holds the sub-class (0 = not expanded), `cluster` the parent;
label_truth_readout.py):
- arms: 92 sub-classes with >= 30 particles; per-arm angle std ratio: best 0.46, median 0.70 (klem leaf states: best 0.74;
  boards' coordinate/angle |corr| <= 0.48). Sub-classes of 90-190 particles localise single arms (arm 1, arm 4 mostly).
- tethered: 97 sub-classes; best site occupancy score 0.30 (= occupancy 0.85 or 0.15 at one site), median 0.62 (klem: no
  region states, |corr| <= 0.08).
=> pooling before the decision (the user's 2D expansion) is what makes a 12 kDa moving domain at this SNR visible; the
per-image 3D E-step is at the noise floor there. Regime split: klem tree + zoom for discrete states of large data,
2D route for small domains / small data (RAG); both can be chained (sub-class labels -> solve3D_cavgs -> states).

**Why:** the first positive on the hard regime, and it comes from the user's own route; it also answers the RAG question.
**How to apply:** for small domains or small datasets start with cls_expansion sub-classes; judge by the per-arm /
per-site readouts; next = solve3D_cavgs on these sub-classes and the arm occupancy of the resulting states.

**Bag prior inside the board (2026-10-17 evening):** the Python board with the cls_expansion sub-classes as bags (`--bags/--bag_parents --bag_beta 2`, parent-centred, arm-focused 3x3 lattice on the 13.8k bagged particles of the light arms set) gives the same per-cell arm std ratio as the per-particle E-step (median 0.85 vs 0.88, best 0.71 vs 0.66) while the sub-classes themselves reach 0.46: sharing the particle log-likelihoods inside a bag does not recover what the averages see, because the per-particle 3D likelihood is the blob. Next variant = the averages-level sub-class E-step (`--avgs --avg_mode pool|unmix`, bag_arms_avgs.sbatch). Scripts: rnasep_klem/bag_arms on batch2.

**Averages-level sub-class E-step (same evening):** `--avgs cls_expansion_cavgs.mrc --avg_mode pool` (mean max r 0.93, confident sub-class -> cell mapping) gives median 0.85 / best 0.72, `unmix` 0.81 / 0.66: no better than per particle. The cells never differentiate on the 12 kDa arm from 13.8k particles, so the confident mapping is view/noise, not arm angle. Transfer has to start from the averages' own latent (cluster the sub-class averages into states, then particles), not from a likelihood against undifferentiated cells. Dump trap: print_project_field on the expansion project needs the kernel snapshot build (~/SIMPLE/build prints nothing, backtrace).
