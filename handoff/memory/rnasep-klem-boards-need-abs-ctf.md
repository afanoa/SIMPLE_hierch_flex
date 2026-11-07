---
name: rnasep-klem-boards-need-abs-ctf
description: TRAP (2026-10-16): the 692k RNase P stacks are phase-flipped, so klem boards must apply |CTF| (fboard CTFABS=1; the Fortran key ctfabs 1); the signed CTF (ctfabs 0, picked by the 16-12 A calibration) made the lattice axis a defocus axis (rho -0.61), a 17 % low-defocus junk streak and a fake low-defocus "family 2" ("three-pronged" maps)
metadata:
  type: project
---

Symptom: parent-board z1 vs defocus Spearman -0.61, junk mass vs defocus -0.76 (17 % of particles junk, mean defocus
0.84 um), families split at 1.24 vs 1.51 um, family 2 matched 0-4 of the 55 deposited apo maps, refined states 9-16 A.
Cause: at 1 um the first CTF zero is at 14 A, inside the 20-10 A bands; with phase-flipped images and a signed CTF
model the sign disagrees past the first zero, more so at low defocus. The 16-12 A correlation calibration that chose
ctfabs 0 was not sensitive enough (SUPERSEDES that line in [[rnasep-klem-real-data-setup]]).
Fix (f692ca / f692caps, 18:52): ctfabs 1 -> rho -0.11, junk 8 %, families at equal defocus; + pscale 1 -> rho -0.05,
junk 4 %, parent cells identify 9/10 (0.41) and 18/55 (0.49) vs 5/10 (0.33) before.

**How to apply:** every real-data board: fboard.sh CTFABS=1 PSCALE=1 (chain8 launched 18:55 as f692c with these);
the f692 landscape (47/62/104/180) and f692b are VOID; toys (simulated with the full CTF) stay at ctfabs 0. Check the
stk segment's ctf= key before trusting any CTF convention. Related: [[rnasep-denoised-stacks-are-phase-flipped]],
[[rnasep-f692-landscape-and-family2-defocus-2026-10-16]].
