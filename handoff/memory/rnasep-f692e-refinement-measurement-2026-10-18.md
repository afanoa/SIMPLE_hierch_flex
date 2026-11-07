---
name: rnasep-f692e-refinement-measurement-2026-10-18
description: how far refine3D_states takes f692e states (4 states, 38k/15k/5k/4.9k-closed, hard labels, local 20->5 A from consensus poses, 26 min): unfiltered half-map FSC 5.3->5.4, 10.2->10.2 (junk-suspect), 7.8->7.8, closed 7.6->6.0 A; closed occupancy 0.94->0.84
metadata:
  type: project
---

refine_states_f692e.sbatch (job 53958452, 2026-10-18 10:52-11:17, snapshot 5e1e6aeee, 20 parts x 8 thr, 8 its):
unfiltered half-map FSC before -> after (SIMPLE's filtered halves in brackets): state 1 (37.9k, leaf f1__1_s01) 5.3 ->
5.4 [5.1]; state 2 (15.3k, f2__1_s08, junk suspect: no deposited match 0.155) 10.2 -> 10.2; state 3 (4.9k, motion_1_s01)
7.8 -> 7.8 [7.6]; closed zoom state (4.9k) 7.6 -> 6.0 [5.7], closed-arm occupancy 0.94 -> 0.84 (consensus-started
local search dilutes the arm). Reading: the leaf states' poses from the consensus refinement are already at the
particle-count limit for the large states (no gain); a small PURE state gains ~1.6 A; a junk-like state gains nothing;
hard labels only (the soft-weights build is not in this snapshot). User: "I wouldn't pursue refinement just yet" then
"try 3-4 states with varying populations". Related: [[rnasep-f692-corrected-landscapes-2026-10-17]].

**Ceiling of the data as extracted (2026-10-18):** the lab's 25_refine3D_states on the 692k set (box 160 @ 1.464 A, no
CTF refinement / polishing) reaches 4.4-4.8 A for its largest states (FSC 0.59-0.61 at 6 A, 0.25 at 5 A, 0.06 at 4 A).
Scaled by the paper's B (170 A^2) from ~300k to 38k particles that predicts ~6.3 A; our 38k state at 5.3 A is AT/above
the ceiling of this extraction. Paper-level (2.9-3.9 A) needs the pre-processing (0.73 A/px box 480, CTF refinement,
polishing) on the pure subsets, not more from the heterogeneity method; the landscape is decided at 8 A anyway.

**Rigid-core alignment attempt (core_align.sbatch, 2026-10-18 17:xx):** this snapshot's refine3D has no mskfile, so
the substitute was refine3D (refine=neigh, lp 6, 6 its) against the family parent multiplied by a core mask (envelope
minus the moving regions, 1.1 % of the box). Result: the alignment got WORSE (family-1 largest state global FSC 5.3 ->
6.3 A, arm-region 3.8 -> 5.3 A; the closed sub-state's occupancy readout collapsed to 0.14, both maps weak): a masked
reference with an unmasked particle comparison treats the periphery as noise and drifts the poses. Proper focused
alignment needs real mask support in the refinement (newer builds; record-width trap across builds). Also the test
picked two family-2 states that are absent from family 1's project (only states 1 and 4 reconstructed). Parked.

**Rigid-core alignment WORKS with pcg_mskfile (user's pointer; core_align.sbatch 2026-10-18 18:xx):** refine3D of all
family-1 particles against the family parent with pcg_mskfile = core mask (envelope minus moving regions; the PCG
support projection solves the reference inside the core every iteration), refine=neigh lp 6, 4 iterations, 20 x 8,
~30 min. States re-reconstructed at the core-aligned poses (no per-state refinement): global FSC 5.3 -> 5.0 (38k),
7.6 -> 5.9 (8.6k), 7.8 -> 7.1 (5.3k), closed 7.6 -> 6.2 A with occupancy 0.92 kept; arm-region FSC 5.7 -> 5.2, 6.7 ->
6.3, 4.3 -> 4.0 (largest state 3.8 -> 4.7). Trap: the first masked-REFERENCE substitute (no pcg_mskfile) made poses
worse; the family project must be made single-state (selection all 1) before refine3D (10-state docs crash the parts);
my ITERATION/RESOLUTION log checks read stale logs (execution dir numbering continues: 3_, 5_refine3D). Next: family 2,
then rebuild the 63 leaves at core-aligned poses (consolidate with the core-aligned projects).
Family 2 core-aligned too (5 its, 23 min). The 63 leaves rebuilt at core-aligned poses (pooled_dd TAG=core_region,
PS_F1/PS_F2 overrides in consolidate.sbatch): FSC median 9.4 -> 7.3 A, best 5.3 -> 4.8, worst 11.2 -> 10.2, states at
<= 7 A 10 -> 23; docked 10/10; strict recovery 3/55 (was 7/55; mean best-match corr 0.36 vs 0.39: the core-aligned
frame matches the docked conformer deviations slightly less - check frame/periphery before reading it as a loss).
Delivered as ~/Downloads/.../f692e_leaf_states_corealigned/. Rigid-core alignment = one refinement per family, the
cheapest resolution gain found; it belongs in the chain after the family refinement.
Frame drift confirmed and fixed (core_frame_check / core_dock_states.sbatch): family-1 map at core-aligned poses vs
toy consensus 0.930 (family poses 0.983); dock_volpair of every core-aligned state onto the consensus -> 0.984; strict
recovery on the DOCKED core-aligned set 8/55 (family poses 7, undocked 3), mean best-match corr 0.44 (0.39 / 0.36),
docked 10/10. Rule: after any whole-set alignment, dock the maps (or transform the poses) back into the docked frame
before deposited-conformer readouts. Delivered folder f692e_leaf_states_corealigned now holds the docked maps.
