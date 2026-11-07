---
name: rnasep-arm-path-maps-two-state-test-2026-10-21
description: "RNase P arm (f692st): two-state rank test on ordered maps; whole-set junk-free (536k) map-level classes closed 9.2 % / intermediate 19 % / open 72 %, rank test fails 7-11x noise (not two-state); WITHDRAWN: r1 lockstep (built into linear-coordinate windows) and the third population (those leaves are junk)"
metadata:
  node_type: memory
  type: project
  originSessionId: f0f8499a-d0a6-450b-9f93-483f98f27006
  modified: 2026-10-22T00:55:27.994Z
---

2026-10-21 arm-landscape task (see [[arm-landscape-deconvolution-truth-tests-2026-10-21]]); maps `view_rflex/landscape_f692st_r1` (lp 10
unsharpened, halves/), runs `rnasep_klem/landscape/real_f692st_r1` (coordinate, rec/, rank.txt).
- **Two-state test** (`land_rank.py`): with two states, every arm sub-region is a linear function of the closed fraction, so the
  (sub-regions x ordered maps) matrix is rank 1. Truth check: two-state mixtures 0.004, rigid swing 0.148. On data the second singular
  value is compared with the (even - odd)/2 noise.
- **r1 windows: closed site in lockstep.** The closed-site shells at 60-118 A decline proportionally (ratio 0.75): no tip-first or
  base-first sequence. The hinge-side 40-60 A shell loses only about 20 %.
- **r1 windows: open site not in lockstep.** With the open site added the ratio is 1.30. The open site fills over windows 1-6 and
  saturates at ~0.75, while the closed site keeps emptying (0.64 -> 0.47): particles leave the closed site to a third position, consistent
  with [[f692st-final-deliverables-2026-10-21]] "~36k arm in neither position".
- **Whole-set populations.** The 64 plain leaf half-maps unmixed (NNLS) onto 63 deposited conformers, classed by the conformers'
  closed-site readout (>= 0.8 / 0.45-0.8 / < 0.45), give closed 0.076 / 0.077 (even / odd), intermediate 0.160 / 0.156, open 0.764 / 0.767.
  - Leaf-level readouts: mean closed-site 0.31; only 3.1 % of particles are in leaves >= 0.8 (r1_s5, r51_s5).
  - Class shares depend on the maps (crop-48 held-out vs SIMPLE full maps differ by up to 0.2 per class). Quote the readouts first.
- **Coordinate.** Inside r1 the start is the cross-half leaf-map PCA (PC1 25 %). For the whole set, PC1 is families, so it starts from
  the leaf closed-site occupancy (`land_real.py --init arm_raw.txt`).
- **The rigid swing of the closed arm is not the open state.** Its open-site occupancy is at most 0.71 and open leaves fit it with
  R^2 0.1-0.5. Several deposited conformers (conf19 / 30 / 50) lie on its early, partially swung part.
- **TRAP:** `arm_readout.py` given a DIRECTORY reads only `state_*` / `recvol_state??` files; pass deposited maps as files.
**CORRECTIONS 2026-10-22.**
- The user says the "arm in neither position" leaves are JUNK: r6_s1/s2/s4/s5/s6/s7/s8, r4_s1/s6/s7, r2_s4 (leaves 11, 24, 25, 27, 34-38,
  40, 41; 107,116 particles; relayed as "about 36k").
- Junk-free: classes 0.091-0.093 / 0.187-0.191 / 0.718-0.720, mean closed-site 0.345, rank test 11.4x / 7.25x.
- The r1 "lockstep" is largely built into windows of one linear coordinate. Windows along an OT morph look tip-first by construction.
  Neither is mechanism evidence.
**How to apply:** answer "is open turning into closed" with the two-state test on ordered maps plus half-map noise, not with conformer
class weights, because NNLS pulls two-state wells inward and so fakes intermediates.
