---
name: klem-solvent-mask-and-supportk-verdicts-2026-10-17
description: toy A/B (dev preset, 150k of the 200k apo/holo toy): the family solvent mask from the feedback cells (derive_solvent.py, 2 % of the box) killed every region (no supported states); supported-K refinement gives K=1-2 which refine3D_states refuses (needs >= 3, and 3 segfaults) -> K = max(supported, 4)
metadata:
  type: project
---

2026-10-17 16:55, toyah_solv0 (no mask) vs toyah_solv1 (SOLVENT=1), chain9 DEV=1 LOCALDEV SUPPORTK=1:
- solv1: solvent mask = union of the 49 feedback cells at 0.15 of max, dilated 10 A: only 2.0-2.2 % of the box (7 % soft) -
  far too tight; all 6 region boards reported "no states". Verdict: the mask as derived removes the molecule, not the solvent;
  if revisited use a much looser threshold (<= 0.05) or an envelope from the consensus, and test on the toy first. SOLVENT=0 stays.
- solv0: 6 states at 0.7, 3 of 10 conformers identified (cell count) - much weaker than the full toy chain (22/48/90/173,
  6-7 found) and the real dev run (8/10). Cause: SUPPORTK gave K=2 (fam1, refinement 38 min for 2 states) and K=1 (fam2,
  refused: "requires nstates >= 3 for flex_pca initialization"). The refinement needs several state-specific references to
  help; f692_family.sh now uses k-means K=4 whenever fewer than 4 states are supported.
- Local-worker DEV loop on the toy: parent 167-186 s, child ~4 min, region boards 12-13 min each sequential at 60k particles.

**Why:** two ideas tested and parked with numbers; the refinement's K is constrained by SIMPLE (>= 4 in practice).
**How to apply:** do not enable SOLVENT=1 as is; keep K >= 4 for refine3D_states; judge toy A/Bs by the merge-tree /
truth readouts, not the cell count.
