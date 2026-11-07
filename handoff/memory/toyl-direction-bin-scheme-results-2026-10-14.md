---
name: toyl-direction-bin-scheme-results-2026-10-14
description: toyL2/L3 (RNase P core + two arms +-30 deg, 40k particles, SNR 0.02/0.05): both latents see the swings (3D per-particle R2 0.93-0.97), direction bins beat 2D classes down to 200/bin, averages-flex stitching weak (pooled particle latent 0.98), refining against MIXTURE references (windows along one axis, mean volume, averages) worsens poses in every mode, but against a STATE-PURE reference poses drop 4.2 -> 1.4 deg: state granularity decides
metadata:
  type: project
---

Toy (2026-10-14, `/mnt/beegfs/elmlund/afan/toyL2`, scripts `scripts/cls_expansion_proto/subclass_pool/toyL/`,
notes section 20): RNase P consensus with region 1 (42 % mass) fixed and regions 2-3 swung +-30 deg (5 levels each,
125 states, far ends 38-40 A), box 96 @ 2.44 A, 320 particles/state, SNR 0.02, CTF; consensus poses = refine3D
neigh from 10 deg off against the true-pose consensus: 4.5 deg median error, state-dependent 3-6 deg.
- Per-particle 3D flex latent (lp 10, 16 comps) orders both arms: ridge R2 0.93 / 0.94, 15-NN level purity 0.83.
  Pooled per sub-class R2 0.98 / 0.99. 2D sub-classes (solve2D 100 classes + cls_expansion pose mean) follow the
  angles: within-class ANOVA F 19 / 46 vs 0.8 permuted. Neither latent is blind to large motions here.
- Direction bins (spherical k-means of refined directions, class = bin, make_cavgs with `nparts=1`, cls_expansion
  unchanged) order the angles better than image classes: F 48 / 140 at 100 bins (400 members), still 21 / 53 at
  200 bins (200 members). Per-bin signal was not the binding constraint.
- Stitching by flex on the sub-class averages (1200 images, 8 comps) gives R2 0.17 / 0.67 only; the pooled
  per-particle latent is the better stitching coordinate on the toy.
- 3D window maps at the consensus poses recover the true window deviation at corr 0.5-0.6 (recovered fraction
  0.55-0.78, amplitude 0.9-1.8), the toy analogue of RNase P's partially reproducible window deviations.
- Per-particle own-map refinement of window halves (refine3D neigh, 5 its, lp 8 or 15; halves of 1k-11k) makes the
  poses WORSE (4.5 -> 5.4-7.9 deg) and does not raise the correlation with the truth (0.62 -> 0.55), even from the
  EXACT reference (-> 5.3 deg); what rises is the reference's imprint (recovered fraction 0.55 -> 0.73) and the
  even/odd reproducibility (0.83 -> 0.92).

- The oracle that explains it (toyL3): the particles of 12 single true states refined against their OWN true volume
  go from 4.2 deg (1.7 central, 5.5-5.8 extreme states) to 1.4 deg (`cont`) / 1.6 deg (`neigh`). The consensus error
  is state bias, and a state-pure reference removes it; mixtures along one axis (thirds, eighths, the mean volume,
  averages blurred over a bin) do not.

- But no reference the data can build helps: own reconstructions of structure-pure groups (320-2700 particles,
  true or latent-assigned) give 4.0 -> 3.8 deg at best and converge there in two passes; docking them onto the
  consensus first makes it worse (-> 5.6 deg). The own map is the state blurred toward the consensus by its members'
  biased poses, so the loop sits at a biased fixed point.

**Why:** refinement against a reference reinforces that reference's deviation; references built from biased poses
(or mixtures, or the mean) carry the bias, so the poses stay while reproducibility and "recovered amplitude" rise.
Only pose error vs truth or correlation with truth certify a gain; a state-pure reference in the TRUE frame is what
moves poses (oracle 4.2 -> 1.4 deg).

**How to apply:** treat the RNase P gold-standard gain (section 19.3) as unproven fidelity. Per-window / per-state
refine3D from own maps is self-consistent and does not move the motion into 3D; the next increments are an explicit
deviation model per state (consensus + modes) or a joint pose-state search from coarse references, both testable
on toyL3 with the toyL scorers (`run_toyL_stateloop.sbatch`, `toy_map_score.py`). refine3D `cont` is the least
damaging mode; `neigh`/`greedy` scatter by degrees; flex_pca's state labels must be cleared by a selection first. Traps: flex_pca leaves state labels in
ptcl3D (run a selection with an all-ones keep file before refine3D); `solve2D` not abinitio2D in the coupled
snapshot; register refined maps before any map comparison. Related: [[rnasep-per-window-poses-gold-standard-2026-10-14]],
[[per-window-refine3d-shared-reference-inflates-reproducibility]], [[toy2-hinge-swing-test-system]].
