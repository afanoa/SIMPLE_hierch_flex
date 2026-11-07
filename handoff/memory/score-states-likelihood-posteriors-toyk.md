---
name: score-states-likelihood-posteriors-toyk
description: "2026-09-28 score_states (private exec, branch flex-weights-refine) gives sharp, near-calibrated image-likelihood state posteriors on toyK (95% vs truth); the chi-square calibration scale works; the prob-search exp(-dist) is on the wrong scale"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-09-28T15:46:24.898Z
---

`simple_private_exec prg=score_states projfile=... nstates=K vol1..volK objfun=euclid sigma_est=global lp= trs= center=no`
(commander in `simple_commanders_prob.f90`, shared memory only): per particle and state, fixed current pose + shift
search, loss re-evaluated with `gen_raw_euclid_vals`, R = L * wsqsum_ptcl_now; scale c = M_eff / median best-state R
with M_eff = pi (k2^2-k1^2)/2 * mask area fraction; r ~ pi exp(-c R), pi by EM; writes flex_weights_state_NNN.bin
(registered) + score_states.txt.

toyK (`/mnt/beegfs/elmlund/afan/rnasep_regions/toyK`, script `~/toyK.sbatch`, ~4 min end to end): toy2 states
1/7/21/27, 5k each, box 96 @2.44, SNR 0.1, poses after 1 greedy iteration from truth, exact state maps as refs.
Result: argmax accuracy 0.951 vs truth, median top posterior 1.00; calibration mildly OVERconfident (bins
0.54/0.47, 0.71/0.63, 0.89/0.81, 0.997/0.980 = mean posterior / accuracy). So image-likelihood posteriors are sharp
where the latent classifier was flat (RNase P median 0.33).

**Scale trap:** the probabilistic search samples states with exp(-(d-dmin)) where d is the particle-power
NORMALISED loss (~1); the code comment calls it calibrated but differences are tiny on that scale, so those
probabilities are near-flat. Any weight derived from prob tables needs R = d * wsqsum and the chi-square c.

Same day: the flex=no refine3D_states control did NOT reproduce the RNase P reference (label ARI 0.045, deviation
corr <= 0.33, final labels weakly separable in its own final-pose latent, Mahal^2 median 0.40 vs 3.74); the user
judged flex=no as not working as expected, so it is not a discriminating control. Reference states 3 and 5 match
no control map (corr 0.79-0.85): the known junk sinks. Related: [[rnasep-final-poses-make-states-separable]].
