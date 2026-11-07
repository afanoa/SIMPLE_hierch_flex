---
name: refine3d-states-nstates3-segfault
description: TRAP: refine3D_states (kernel snapshot 5e1e6aeee) segfaults deterministically with nstates=3 right after the initial state reconstruction (seen twice, local parts); nstates 4-12 fine; also needs maxits >= 5
metadata:
  type: project
---

2026-10-17: f692smoke3 and f692smoke4 family 2 (25k particles, 3 supported states): `refine3D_states nstates=3 pose_policy=local
lpstart=10 lpstop=8 maxits=5 nparts=6` dies with SIGSEGV immediately after `RECONSTRUCTION MASTER PHASE (gridding)` of the
initialisation, both times at the same point; family 1 with 4 or 5 states in the same runs completes. nstates=12 (purify)
and 10 (production) run. Not reported to the user's main tree yet; the chain falls back to consensus poses when it happens.

**Why:** a deterministic crash in the user's SIMPLE build, not in the klem port; it silently costs the family its refinement.
**How to apply:** when the supported K is 3, expect the crash; report it to the user as a SIMPLE bug with the smoke3/4 fam2
refine.log; do not pad K to dodge it.
