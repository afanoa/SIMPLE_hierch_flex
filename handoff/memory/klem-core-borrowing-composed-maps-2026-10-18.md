---
name: klem-core-borrowing-composed-maps-2026-10-18
description: user idea "borrow signal from the parent node": composed state maps = per-family parent (all family particles at the family's refined poses, FSC 4.3/5.0 A) outside the state's region + the state inside (compose_state_maps.py, parent_fam_compose.sbatch); the consensus-frame parent mismatched family 2; FSC of a composed/weighted map is not a resolution (report parent FSC + native in-region FSC)
metadata:
  type: project
---

2026-10-18 13:00-14:00. compose_state_maps.py RUN CONS PARENT{F} OUT: M_h = P_h + m (S_h - P_h) per half, S
scale-matched to P on the support, m = the leaf's region mask (regions{,_motion}/region_N.mrc) low-passed 15 A.
Traps: (1) the lab's rec_final_state01_{even,odd} (25_refine3D_states) correlate to Nyquist -> not independent halves
(composed FSC "3.0 A" everywhere); (2) a consensus-frame parent (all 692k at consensus poses, 4.5 A) mismatches
family 2's refined frame (its cores got worse, 12 -> 15.6 A); per-family parents at the families' refined poses fix it
(fam1 4.3 A, fam2 5.0 A). Result: family-1 leaf cores 6.9-7.6 -> 6.0-6.2 A, small/zoom states' cores 6.3-9 -> 4.3 A,
native in-region FSC 3.9-6 A unchanged (as it must be). User caveat: FSC does not apply well to weighted/composed
reconstructions -> readout = parent FSC (core) + state's native in-region FSC only; soft-weight refinement must use
weights for half A from half-B references. Delivered: ~/Downloads/.../f692e_composed_maps/ (63 leaves + 4 closed
sub-states + fam1/fam2 parents, fsc_table.txt). Related: [[rnasep-f692e-refinement-measurement-2026-10-18]].
