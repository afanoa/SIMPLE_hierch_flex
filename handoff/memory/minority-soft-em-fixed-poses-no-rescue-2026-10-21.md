---
name: minority-soft-em-fixed-poses-no-rescue-2026-10-21
description: "soft (responsibility-weighted) vs hard EM at fixed poses does NOT rescue the RNase P closed minority on real r22; soft end components balloon into sinks (2k -> 40k), hard drains states; klem got the closed state by OVER-SPLITTING (K=30 in a 240k family, ~8k states) + the coordinate pass"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-21T18:44:17.621Z
---

minor_em.py (routes/toy55; runs rnasep_klem/minor_em_r22, minor_em_s55): cross-half references, arm-mask deviations, per-shell
data scale, fixed parent poses, crop 48. Real f692rf r22 (381k, ~2.6 % closed), closed-arm of the prototype's own maps (they read
~0.15-0.35 LOWER than SIMPLE reconstructions of the same labels; compare only within the prototype):
- equal-width arm seeds (end bins 2.1-2.8k/half, clean gradient 0.45 .. -0.05): soft -> end state 40k/half at 0.41 (best 0.42);
  hard (classification EM) -> 5 states drain, best 0.46 at 28k.
- refined (evened) 12 states: soft 0.47 -> 0.55 at 22.7k/half.
- synthetic sim55_200k with noisy-truth seeds (end bin 75 % closed): both scramble in iteration 1 (middle bins identical maps);
  near-consensus components become sinks (35-46k at 15 % closed). Seeds must follow a real data axis for such a test.
Why klem worked (f692r): refine3D_states K=30 inside the 240k family holding the closed particles -> ~8k states, the closed pair
16/15 (5.5k/6.1k, 0.87/0.83) owned states; the 1.03 cell = coordinate pass (all particles scored vs the 30 maps, PCA, corner cell 2.5k).
Hard refinement drifts to ~N/K per state, so a minority m needs N/K ~ m: rflex's K <= 8 and KFLOOR 15000 block it at every node.
Tracer f692rf: closed share 1.8 % (root) -> 2.0 (r2) -> 2.6 (r22) -> 8.2 % (best 57k leaf); synthetic dtree_55_200k_lat best leaf 20 %.
**Why:** the user asked for a better way to get minority states; soft assignment was the first hypothesis and is ruled out.
**How to apply:** pursue depth with arm-focused (cheap-children) masks + small states (KFLOOR 5000; f692st on rflex10, the unit-7
test rflex_f692u7) and over-splitting per node (K ~ N/8000, vol up to 99 with snap focus3); see [[klem-integrated-hierarchical-model-vision]],
[[junk-free-component-needs-scale-fair-maps-2026-10-21]].

SITE ROUTE WORKS (site_r22.sbatch, minor_em3.py --post; rnasep_klem/site_r22, maps view_rflex/site_r22): 2 classes at fixed poses
that differ only inside a closed-SITE mask (unit 7's cheap-children variation mask: closed site 0.68, open 0.01), all 381k r22 particles,
cross-half refs: the class keeps growing (2k -> 65k/half, posterior a continuum, no clean mode) but the RANKING is informative: SIMPLE maps
of the top N by posterior: 2.5k 0.96 / open 0.80, 5k 0.91, 10k 0.86, 20k 0.81 (seed bin 0.55). Open site stays occupied -> closed-site and
open-site occupancy are separate binary features (klem node 54: 4 corners); the proper closed state = closed occupied AND open empty.
Depth route (f692u7, one level below r22 state 2, cheap-children mask on the closed site, K=8 ~7k states): state 6 = 8.8k at 0.87 / 0.69,
but splitbench15's beyond-EM re-test (16.7-10.2 A) merged all 8 states (no power for low-resolution differences) -> rflex11 DECIDE=all.
BEST SO FAR (17:40): f692st (rflex10 root with the broad flex-state mask, then rflex13: falloff masks, KFLOOR 5000, splitbench16 --focus,
no borrowing) node r1 (depth 1 = root state 5, 92.6k, 37.6 % of the old closed) -> refined state 5: 14,775 particles, closed-arm 1.00 /
open-arm 0.47 (= klem's corner 1.03/0.55 at 6x the particles), 2,843 old closed (43.8 % of the node's). The proper closed state at DEPTH 1.
