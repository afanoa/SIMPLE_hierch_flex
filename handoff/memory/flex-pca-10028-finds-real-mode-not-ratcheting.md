---
name: flex-pca-10028-finds-real-mode-not-ratcheting
description: "With DEFLATE_BG on, flex_pca finds a highly reproducible intramolecular mode on EMPIAR-10028 (cross-fit cosine 0.89-0.98) but rigid-body ratcheting is NOT measurably present in the latent"
metadata:
  type: project
---

Campaign of 20+ runs on EMPIAR-10028 (P. falciparum 80S, 83626 ptcls, 2026-08-21).

**The positive result.** With `SIMPLE_COV_DEFLATE_BG=1` the leading component is a
genuine mode, not solvent:

| component 1 | baseline | after deflation |
|---|---|---|
| power INSIDE molecule support | 0.234 | **0.767** |
| power in top 10% of voxels | 0.573 | 0.889 |
| corr with consensus (scale) | 0.221 | 0.002 |
| signedness | 0.000 | 0.000 (dipolar = mass transport) |

Six independent fits differing in neigs (10/20/40), lp (16/20 A) and mskdiam
(300/320) recover the SAME direction at **|cos| 0.89-0.98**, against a null of 0.042
(PC1 vs high-index components). Baseline PC1 correlates only **0.10-0.14** with all of
them -- deflation did not improve PC1, it REPLACED it. Not a knob artifact.

**The negative result.** Rigid-body ratcheting is not detectable in the per-particle
latent at this operating point. `z` of the rotation direction vs random latent
directions, across every configuration tried:

    baseline +0.11 | dfl4bg +0.73 | n10_dfl4bg -0.51
    p2_lp16_n10 -0.39 | p2_lp16_n40 +0.99 | p2_lp16_n5 -1.06 | p2_lp18_n10 -0.57

Pure scatter around zero, while the BASIS improved ~4x (`R2_k5` 0.068 -> 0.287,
nuisance ratio 5.36 -> 2.05, eigenvalue-order violations 15 -> 1). Representation
improved; the per-particle measurement did not. Ruled out along the way: band
(82% of the rotation template survives lp=20), rank
([[flex-pca-transport-rank-is-not-the-barrier]]), placement (GMM=0 identical),
and detector shape (subunit-local templates score LOWER than whole-volume).

Operating point vs cryoSPARC 3DVA on the SAME dataset, which does report ratcheting:
3 A consensus + full-res particles + ~5 A filter, vs our 4.05 A consensus,
box_crop=64 (7.54 A pixels) and a 20 A band. 3DVA's tutorial also uses a curated
subset where ours is all 83626 particles from the abinitio3D.

Instruments on fwl-c142339: `score.py` (per-run table), `latent_rot.py` (the decisive
z), `span.py`, `nwidth.py`, `lobe_rot.py`, `what_is_pc1.py`, `pc1_repro.py`,
`mask_tradeoff.py`, `defl_safety.py`. Results in `~/campaign_results.txt`.
See [[flex-pca-nuisance-deflation-flags-unset]] for the flag attribution.

## 2026-08-22: CLEAN RERUN -- the ratcheting negative now stands, the positive does NOT

Everything above was measured on a particle set built by a BROKEN class selection
(see [[map-cavgs-selection-ranked-stack-trap]]): classes were chosen by INDEX, so
rank-1 was discarded and 46 of the bottom 60 kept, including classes at 68-121 A.

Rebuilt from the pre-selection project snapshot (all 105,247 particles ON -- note
`map_cavgs_selection` only turns particles OFF, never back on, so re-selecting on an
already-selected project leaves particles permanently stranded), correct rank-based
cut, then a full from-scratch `abinitio3D mskdiam=320 pgrp=c1`:

    100,667 particles | consensus 4.05 A (even/odd FSC) | mean-scale OK on all configs

| config | ratio | R2_k5 | rot_state med | latent z |
|---|---|---|---|---|
| n20_lp16 | 2.44 | 0.081 | 0.033 | -0.33 |
| n10_lp16 | 2.31 | 0.055 | 0.041 | -0.79 |
| n40_lp20 | 2.83 | 0.039 | 0.020 | +0.47 |
| n20_lp20 | 2.47 | 0.054 | 0.028 | -0.50 |

**The ratcheting negative is now TRUSTWORTHY** -- z spans -0.79..+0.47 on a correct
particle set with a correct consensus. Every prior null was on contaminated data.

**The positive above is SUSPECT.** `R2_k5` fell from 0.287 (contaminated) to 0.039-0.081
(clean) at matched settings, and `rot_state` med from 0.093-0.122 to 0.020-0.041. The
earlier "mode reproducible at |cos| 0.89-0.98" was very likely driven by the junk classes
and has NOT been re-verified on clean data. Do not cite it until it is.

Component 1 does remain genuinely above the noise floor (ratio 2.3-2.8, eigviol down to 1
at neigs=10), so flex_pca is fitting real structure -- just not rotation.

Eliminated as explanations: band, rank, placement, detector shape, nuisance deflation,
mask, particle selection, consensus quality. The one untested difference from cryoSPARC
3DVA (which reports ratcheting on this same EMPIAR entry) is FIT RESOLUTION: we fit at
box_crop=64 -> 7.54 A pixels; 3DVA works at full resolution. Verified project for that
test: `/usr/local/data/atifao/10028_clean/1_selection/1_abinitio3D/`.

## 2026-08-22 later: the clean negative is PROVISIONAL again -- n_probe_iters confound

prod10 (the original production run) used `n_probe_iters=25`. The A/B harness cut it
to 3 for cycle speed and the cut silently propagated into EVERY subsequent harness:
campaign, phase 2b, the overnight clean run, bc128. So all post-prod10 results ran at
3 EM iterations, not 25 -- under-convergence is a live confound for the ratcheting
negative (a weak mode separates from the noise bulk slowly in power iteration).
Re-testing at 25 in `~/run25.sh` on the box (probe25 at bc64, then bc128 n10/n20).
Lesson: when deriving a fast harness from a production script, diff the FULL argument
list, not just the env block.

Also staged: `github.com/zhonge/cryodrgn_empiar` cloned to
`/usr/local/data/atifao/cryodrgn_empiar/` -- cryoSPARC-refinement poses for ALL
105,247 particles of this stack, parsed at D=360 (our box), stack order = import
order. Converting these to a SIMPLE project (Euler + PIXEL shifts, gate on even/odd
FSC reaching ~3.x A) and running flex on them is the decisive H1-vs-H2 split: rotation
appears -> our poses absorb it; still null -> our estimator.

## 2026-08-22 RESOLUTION: the negative was UNDER-CONVERGENCE. n_probe_iters=25 is load-bearing.

Clean A/B, identical config (10028_clean, n20/lp16, deflation 4+BG, polar), only the
EM iteration count differs:

| | R2_k5 | R2_k10 | R2_all | rot_state max/med | mono | latent z |
|---|---|---|---|---|---|---|
| 3 iters (the whole campaign) | 0.081 | 0.216 | 0.269 | 0.098 / 0.033 | 0.077 | -0.33 |
| **25 iters (prod10 setting)** | **0.300** | **0.493** | **0.672** | **0.478 / 0.315** | **0.385** | +0.89 |

Rotation content ~10x at basis and state level; delivered states now correlate with
rigid-rotation templates at 0.32 median / 0.48 max with a monotone state progression.
The ratcheting "negative" of 2026-08-21/22 is RETRACTED: every one of those runs used
n_probe_iters=3 (the A/B harness cut 25 -> 3 for speed and it propagated silently).
A weak mode needs many power iterations to separate from the noise bulk; 3 was starving
it. H2 (estimator under-convergence) was the answer, caught by the user remembering
"we had default 25 probe iters".

Caveat on the latent z (+0.89, not yet >2): the z statistic is INTERNAL-relative --
its random-direction null lives in the same latent span, so as the basis aligns with
rotation the null inflates too. Use rot_state/R2 (external templates) as primary once
the basis is converged; z is most useful when the basis is rotation-blind.

Still in flight to consolidate: bc128/lp10 at 25 iters, and the cryoSPARC-pose gate
(pose-quality comparison, still worth having). Rules going forward:
- **never run flex_pca A/Bs below n_probe_iters=25**; iteration count is a science
  parameter, not a speed knob. If a fast harness must cut it, it measures SPEED only.
- deflation at 25 iters decays to ~2.3%/iter (healthy), vs the 5-6% rank-1 plateau.

## 2026-08-22 FINAL VERDICT (tail-reconstruction experiment): the E-step per-particle
coordinate does NOT encode the rotation -- genuine method gap upstream of delivery

On the VALIDATED cryoSPARC-pose project (their refinement poses reconstruct to 3.42 A
through our reconstructor; cryoDRGN and 3DVA demonstrably recover ratcheting from these
exact poses):

- flex latent detection: z = +4.56 / +4.70 (only >4 sigma results of the campaign)
- but state sweep collapses to 0.33 deg, no basis component is rotational (R2_k5 0.013)
- **decisive**: sorting particles by our per-particle rotation coordinate s_i and
  directly reconstructing the extreme deciles (10.5k each, 3DVA-"intermediates" style,
  no averaging/shrinkage) gives **hi-vs-lo = 0.52 deg, corr 0.141** -- the tails do not
  rotate. With per-particle SNR ~0.65 along the direction, honest dilution predicts
  ~2-3 deg if s_i encoded the true ratcheting; 0.52 deg rules that out.

So: the significant variance along the template direction is NOT (mostly) the ratcheting,
and the per-particle E-step readout at box_crop=64/lp=16 (~2.5k coefficients onto a
reduced basis) loses information that cryoDRGN extracts at D=256 with an 8-dim latent.
On our own (4.05 A) poses the 4.4-deg state sweep is partly pose-error variance
masquerading as rotation (alignment errors ARE small rotations) -- treat it as an upper
bound, not a capture.

Also refuted en route: "pose absorption ate the amplitude" (a global pose cannot absorb
DIFFERENTIAL rotation); state_axis=1 rendering (leading component is not rotational on
clean poses: sweep 0.14 deg).

Next discriminator (plan doc E1-sharp): a DIRECT per-particle matched filter -- project
t_rot at each particle's pose (full band, CTF-weighted), correlate with the residual
image, no reduced basis. If THAT sorts rotating tails -> the reduced E-step is the loss
and the fix is band/weighting in the E-step. If even the direct matched filter fails ->
single-image SNR for this mode is below per-particle detectability and the readout needs
joint/amortized estimation (what 3DVA and cryoDRGN both do).
Artifacts: ~/tails.sh, /usr/local/data/atifao/10028_tails/{lo,hi}, rot coordinate from
b_cs_n10. Volumes for visual check: tails lo/hi recvol_state01.mrc.

## 2026-08-23: the captured mode is GLOBAL rotation (pose jitter), and DEFLATE_POSE is implemented

User observation from the delivered volumes: our states differ by a WHOLE-MAP rotation
offset, while 3DVA/cryoSPARC rotations are subunit-specific. That identifies the captured
mode as per-particle alignment jitter (a pose error IS a small global rotation), which
coherently explains: the 4.4-deg sweep on our 4.05 A poses (jitter variance), its collapse
to 0.33 deg on 3.42 A cryoSPARC poses (jitter gone), and the failure of tail sorting
(s_i tracks jitter, not conformation). The whole-volume rotation template cannot separate
global from differential rotation -- diagnostics must use the two side by side.

**Implemented `SIMPLE_COV_DEFLATE_POSE=1`** in the box tree
(~/SIMPLE_polarem/src/main/flex/simple_flex_pca_columns.f90, backup
~/columns_backup_posedfl_applied.f90, build rc=0): appends the 6 global pose-derivative
templates (translation gradients + whole-map rotation derivatives about x/y/z, central
differences on the consensus, soft-masked, slots ndfl_sh+1..+6, Gram-Schmidt handles
orthonormalisation). Agnostic and knob-free: pose jitter is nuisance in every dataset
regardless of heterogeneity type; differential (subunit-local) motion is orthogonal to
whole-map rotation and survives. A/B queued (`~/pose_ab.sh` -> poseab_results.txt):
  - our poses: whole-map rotation content in states should COLLAPSE (fake 4.4 deg dies)
  - cryoSPARC poses: does z=+4.6 survive? survives -> genuinely differential; dies -> jitter too.
Standing constraint reminder: templates/tail-sorting are INSTRUMENTS only; shippable
levers are operating point, weighting, iterations, pose/BG deflation, excursion delivery.

## 2026-08-23 CLOSING STATE: DEFLATE_POSE validated; detection is real and differential;
## per-particle fidelity is the one remaining gap

A/B results (poseab_results.txt, csbc128_results.txt on fwl-c142339):

- **DEFLATE_POSE validated.** Our-poses arm: whole-map rotation R2 0.672 -> 0.000,
  rot_state 0.478 -> 0.057, sweep 3.16 -> 0.37 deg. The entire "ratcheting" captured on
  our 4.05 A poses was pose jitter (user identified it visually: states were globally
  rotated consensus copies). Ship `SIMPLE_COV_DEFLATE_POSE=1` as a default candidate
  with DEFLATE_BG.
- **The clean-pose detection is genuinely differential.** cryoSPARC-pose arm with pose
  deflation: z = +4.84 (up from +4.56), and the lobe test flips for the first time:
  subunit-local rotation R2 0.101 vs whole-map 0.001.
- **Band does not buy per-particle fidelity.** bc128 on clean poses: z +3.55 but tail
  reconstruction unchanged (0.52 deg, corr 0.14, same as bc64).

Final gap, precisely: single-pass per-particle projection saturates at SNR ~0.65 along
the mode -- detects (~5 sigma) but cannot sort or render amplitude. Excluded as fixes:
band, iterations (25 needed vs 3, but saturated), delivery rendering, all nuisance
(solvent/envelope/pose). The missing capability vs 3DVA/cryoDRGN is JOINT
eigenvolume+latent estimation at full band with per-coefficient CTF^2/sigma^2 weighting.
Design problem for the GPU path under the real-time constraint.
Instruments (never ship): rotation/lobe templates, tail sorting. See
[[flex-pca-must-be-knob-free]].

## 2026-08-23 (afternoon): weighting certified + band exhausted -> the failure is the M-STEP

Roadmap execution on the 21k dev harness (stride-5 of the cryoSPARC-pose set):

- **M1 audit result: the per-coefficient weighting is EXACTLY the matched filter.**
  `gen_fplane4rec` stores y/sigma(k) and transfer CTF(k)/sigma(k) per coefficient
  (per-particle, colored, phase-flip handled with |CTF|/CTF^2). G and b inherit
  symmetric weights. No weighting bug exists; M1 narrows to the shrinkage/debias export.
- **M2 result: band does NOT fix per-particle sorting.** bc256/lp6 (~20k coefficients,
  3DVA's regime, neigs=8, 25 iters, full deflation): tails 0.48 deg / corr 0.105 --
  IDENTICAL to bc64 and bc128. z=+3.24.
- **The localization:** at the full band the basis spectrum INVERTS (comp1_snr 9.2 vs
  plateau 12.5, ratio 0.74) and per-iteration deflation strips **14%** forever (vs 2-3%
  at bc64) -- the wide band admits envelope nuisance faster than rank-4 shells remove it,
  and rotation R2 in the delivered basis is 0.000-0.035. Sorting fails because s_i
  projects onto a direction that reconstructs ~3% of the rotation. **The E-step is fine;
  the M-step never puts the differential mode into the basis.** 3DVA's remaining
  difference is its clean closed-form PPCA M-step vs our modified EM (per-iteration
  deflation, ONLINE windowing, MCFA mixture, RELPRIOR).
- Dev harness validated: 21k particles reproduce all full-N phenomenology at 1/5 cost
  (~50 min at bc256, ~6 min at bc64); memguard (15 GB floor, enumerated-pid kills) armed.

Next instrument in flight: spike-in benchmark (M0b) -- KNOWN 4-deg lobe-rotation mode
injected into real images (t_diff per-radian template, GT-tail gate self-verifies the
injection). The M-step question becomes directly testable with a ~10-min cycle:
does the EM discover an injected mode (R2(t_diff in basis))? Iterate M-step variants
(deflation rank/schedule, EM_MIX off, RELPRIOR off, plain-PPCA mode, init) against it.
Files: ~/spikegen.py, ~/inject.py, ~/spike.sh; series delivery prototype ~/series.py+sh.

## 2026-08-23 (evening): spike benchmark LIVE -- first floor measured

Spike-in benchmark operational on the 21k dev harness (gate |corr|=0.876 validates the
injection chain end-to-end; python-vs-SIMPLE CTF global sign absorbed into the labels).

**Baseline floor at theta_rms=4 deg (population ~+/-8 deg), standard pipeline (bc64/n20/
lp16/p25, full deflation): R2(t_diff in basis)=0.015, r(s_i,c_gt)=0.501.**

The split finding: the M-step still does NOT discover the injected mode (R2 0.015 --
same blindness as real data), yet per-particle sorting through the near-blind basis
reaches r=0.50 (the 3DVA-level threshold). Hypothesis this supports quantitatively:
the real ratcheting at theta_rms ~1.5-2.5 deg would give r ~0.2-0.3 (r roughly linear
in amplitude) -- tails enriched too weakly to survive reconstruction averaging, which
matches every real-data observation. The 2-deg rung (queued) tests this directly.

In flight: M-step variant sweep at fixed 4 deg (EM_MIX off / RELPRIOR off / ONLINE off /
all-deflation off), each scored R2+r in ~8 min. Instruments: ~/inject.py ~/spike.sh
~/floor_eval.py ~/mstep_sweep.sh ~/rung2.sh; results ~/spike_results.txt,
~/mstep_results.txt. Traps fixed en route: fd exhaustion (memmap cache cap 400 +
ulimit 8192), bare-filename star paths skip rootdir resolution (stack must live in a
subdirectory so the name carries a slash).

## 2026-08-23 (night) BREAKTHROUGH: the DEFLATION eats the differential mode

M-step variant sweep on the 4-deg spike (each run identical except one switch):

| variant | R2(t_diff in basis) | r(s_i, c_gt) |
|---|---|---|
| baseline (all mods on) | 0.015 | 0.501 |
| EM_MIX off | 0.015 | 0.501 |
| RELPRIOR off | 0.015 | 0.465 |
| ONLINE off | 0.015 | 0.501 |
| **ALL DEFLATION off** | **0.458** | **0.932** |

**With deflation off the M-step DISCOVERS the mode and sorting hits r=0.93.** The
orthogonality argument behind DEFLATE_POSE was wrong for off-center lobes: a lobe
rotation about its own centroid CONTAINS a whole-map-rotation component (orthogonality
only holds for symmetric mass), so per-iteration hard projection strips genuine
differential signal every pass and the mode never accumulates. The same overlap logic
applies to consensus-shell deflation (a lobe motion overlaps consensus-shaped
directions). MCFA / RELPRIOR / ONLINE are all exonerated.

**Amplitude hypothesis also confirmed:** r(2 deg) = 0.234 (sign was a label
bookkeeping slip; magnitude is the measurement) vs ~0.25 predicted by linear scaling
from r(4 deg)=0.50. So on real data (theta_rms ~2 deg) the standard pipeline sits at
r~0.23 = unsortable, exactly as observed.

Design implication: deflation must stop being a hard per-iteration projection of the
BASIS. Candidates: deflate only at init; soft prior down-weighting of nuisance
directions; or fit nuisance coefficients JOINTLY with the basis (nuisance as explicit
extra components that absorb their variance without constraining the signal subspace).
The nuisance problem it solved is real (solvent 27x etc.) -- the fix must keep nuisance
OUT of the leading components without projecting signal overlap away.

In flight: (a) attribution -- POSE-only-off at 4 deg; (b) 2-deg rung with deflation off
(projected r~0.45 = sortable); (c) REAL 21k cs-pose data with deflation off.

## 2026-08-23 day-end ledger

Attribution complete: **DEFLATE_POSE alone eats the differential mode** (POSE off,
BG+env on: R2=0.46 |r|=0.94 at 4 deg -- identical to all-deflation-off). BG+envelope
deflation exonerated. At 2 deg (physiological) with deflation off: **R2=0.32, r=0.75 --
sortable**. The benchmark calibrates the sortability threshold: basis R2 ~0.3 needed.

Real data (21k cs poses, BG+env on / POSE off): z=+7.07 (campaign best), lobe1 R2
0.158-0.167 with REPRODUCIBLE axis (+0.80,-0.55,-0.24) across two fits -- but R2 is
half the sortability threshold, and every tail readout (whole-map, unmasked-lobe,
masked-lobe directions; per-component series; varimax axes) stays at ~0.5 deg. Fully
consistent story: true mode ~1.5-2 deg with basis capture ~0.15 -> r~0.3 -> unsortable.
Varimax: all basis axes are already maximally localized compact modes (0.9 top-decile
power); the dominant real heterogeneity is localized density variation (axis 1 = 50%
of subspace variance), rotation lives only as a weak combination.

Overnight run in flight: **bc256/lp6 WITHOUT pose deflation** (m2n) -- the only untested
cell (the previous bc256 had the signal-eater on). If full band + no-POSE lifts basis
R2 toward 0.3, tails become sortable and the series delivers the motion.
Ship-list so far: DEFLATE_BG=1 + EM_DEFLATE=4 stay; DEFLATE_POSE must NOT default on
(eats differential signal; its pose-jitter problem needs the softer joint-nuisance
design); n_probe_iters=25; benchmark harness is the permanent regression instrument.

## 2026-08-24 MAJOR CORRECTION: the template metric was wrong; real localized 40S
## rotation IS being captured

User inspected the delivered volumes and saw the small subunit clearly rotating,
contradicting my lobe-template R2 of 0.124. They were right. Replaced the template
correlation with a TEMPLATE-FREE instrument (`~/local_motion.py`, `~/displacement.py`):
patch cross-correlation -> displacement field -> split patches into MOVING vs STATIC by
|d| -> fit a rigid rotation to the moving set about its own discovered centroid, with
the static set as control.

| run | rot | fit R2 | body centroid | STATIC control | verdict |
|---|---|---|---|---|---|
| probe25 (our poses, jitter present) | 1.18 deg | **0.788** | centred | **0.41 deg, 1.12 A** | GLOBAL = pose jitter |
| DIAG (cs poses, no POSE dfl) | 1.20 deg | 0.246 | **off-centre** (+22,-12,0) | **0.01 deg, 0.29 A** | LOCALIZED = real |
| COUPLED (cs poses, full KxK) | 1.37 deg | 0.268 | off-centre | 0.03 deg, 0.38 A | LOCALIZED = real |

**Key inversion: LOW rigid-fit R2 with a genuinely static remainder = REAL localized
motion; HIGH R2 with a co-rotating "static" set = FAKE global motion (pose jitter).**
The old whole-volume template rewarded global coherence, so it scored the jitter run
highest and the genuine subunit rotation at 0.12. Every "no rotation" verdict derived
from template correlation must be re-checked with `local_motion.py`. The 2-means lobe
split (53/47) was never the real subunit either.

Measured on cs-pose data: ~1.2-1.8 deg subunit rotation, 3-5 A peripheral displacement,
10x contrast vs a static remainder, axis (+0.03,+0.6,-0.8) reproducible across three
independent fits. Short of the literature 6-10 deg -- amplitude compression remains
(state 1 has only 1362 ptcls; 3-state readout is coarse) -- but this is capture, not
absence.

**T0 coupled-solve test: NEGATIVE, ground-truth confirmed.** Full per-voxel KxK
(`COUPLED_DIAG` unset, `rows=210 (full)`, 0.30 GB, no measurable slowdown at K=20/bc64)
vs diagonal, on the 2-deg spike: R2 0.263 vs 0.321, r 0.742 vs 0.753 -- marginally
WORSE. Cross-component coupling is not the limiter; keep `COUPLED_DIAG=1`; the cheap
approximations (rotate-then-diagonal, orthonormal-basis EM, shell-wise KxK, Jacobi)
are moot.
