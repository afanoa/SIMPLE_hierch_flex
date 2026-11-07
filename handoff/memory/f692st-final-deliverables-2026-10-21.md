---
name: f692st-final-deliverables-2026-10-21
description: "f692st state tree (rflex13, real RNase P 692k) finished: 64 leaves / 643k; plain final maps view_rflex/f692st_final_emp2 (+units_index.txt); closed r1_s5 1.07/0.51 with the tip to ~105 A; r51_s5 = a different arm position; 107k 'neither position' leaves = junk per the user (sorted into r6, not flagged by the junk rule); score space = one continuum"
metadata:
  node_type: memory
  type: project
  originSessionId: 25d47be4-7d1d-48d8-abdf-97b1c3c946ab
  modified: 2026-10-22T00:24:37.723Z
---

Run `rnasep_klem/rflex_f692st` (rflex13, falloff masks, focused stage A/B, splitbench16 --focus) finished 2026-10-21 19:39: 64 leaves over 643,140
particles in 11 nodes (r1-r8, r21, r51, r52), no junk leaves. Final job `routes/toy55/rflex_final2.sbatch` (54054943, 38 min):
- PLAIN leaf maps (one reconstruct3D per node, state = leaf, node stage-B poses) -> tw_units (unit = leaf): tree R^2 0.659 < 0.70 -> emp2
  delivered in `view_rflex/f692st_final_emp2` (+ `units_index.txt`: unit -> leaf, particles, arm readouts); raw in `tw_units/f692st_plain/maps_raw`.
  Borrowing keeps the minority states (emp2 vs raw within 0.01: r1_s5 1.058/0.504, r51_s5 0.980/0.791); tree2 moves r5_s5 0.65/0.42 -> 0.57/0.53.
- Closed leaf r1_s5 (14,775): plain 1.07 / 0.51, closed-arm density 0.96/0.87/0.79/0.67 at 90-110 A (the focused map 0.58/0.39/0.23), the
  open-arm region empty at 55-75 A -> a genuinely closed state.
- r51_s5 (4,953): 0.98 / 0.79 but NOT the deposited closed state: outer closed region over-dense (1.1-1.4 at 70-95 A), inner low, open region
  still occupied -> another arm position. ~36k particles (r6_s4/s7/s8, r4_s6) have the arm in NEITHER deposited position (closed ~0.1, open ~0.5-0.6).
- USER (22:45): the 'arm in neither deposited position' leaves (r6_s1/s2/s4/s5/s6/s7/s8, r4_s1/s6/s7, r2_s4) LOOK LIKE JUNK; nice that the
  tree sorted most of them into one family (r6), but the split test's junk rule (relative SSNR < 0.35) did not flag them.
  CORRECTION (motion session): the user's list (leaves 11, 24, 25, 27, 34-38, 40, 41) = 107,116 particles, not ~36k; closed-site 0.135,
  99.6 % open-class. Junk-free set 536k: closed-class 9.2 %, intermediate 19 %, open 72 %.
- Score space (`score_embed_f692st`, figure `viz_f692st/f692st_scorespace.png`): own leaf ranks first 56.9 % (top 3 83 %), PC shares 0.22/0.14/0.11;
  one continuous cloud, closed-arm rising toward one side; r1 and r5 families overlap there; old f692r closed (15+16) concentrate there.
- arm_readout's call column is unreliable at these values (r5_s5 "closed", r1_s5 "intermediate" because open 0.507 > 0.5): report numbers.
Related: [[focus-mask-truncates-mobile-tip-falloff-prior-2026-10-21]], [[minority-soft-em-fixed-poses-no-rescue-2026-10-21]].
**CORRECTION (user, 2026-10-22, relayed by the 2D-to-3D session):** the leaves with the arm "in neither deposited position" are JUNK:
r6_s1/s2/s4/s5/s6/s7/s8, r4_s1/s6/s7 and r2_s4 (107k particles), most of them in family r6. The split test's SSNR < 0.35 junk rule did not flag
them. Exclude them from occupancy and free-energy estimates.
