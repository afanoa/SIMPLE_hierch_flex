---
name: simple-star-shifts-are-pixels
description: "SIMPLE's ori x/y shifts are PIXELS, but the star import reads rlnOriginXAngst verbatim without dividing by smpd — write shifts pre-divided"
metadata:
  type: reference
---

When building a RELION 3.1 star to import into SIMPLE, write the origin columns in **PIXELS**,
even though the column is named `_rlnOriginXAngst`. SIMPLE's internal `x`/`y` are pixels, and
`simple_starproject.f90` maps `rlnOriginXAngst -> x` with NO `mult=` factor (contrast
`rlnPhaseShift`, which does carry one). Export has the same omission. Both directions are only
correct at smpd = 1.0.

Measured on EMPIAR-10049 (108,544 ptcls, smpd=1.23, box 192) — two projects identical except for
shift units, same poses, both through `bootstrap_rec3D`:

| shifts written as | FSC=0.5 | FSC=0.143 |
|---|---|---|
| Angstrom (matching the column name) | 6.95 A | 5.02 A |
| **pixels (pre-divided by smpd)** | **4.82 A** | **4.07 A** |

**Why:** reading the code gives a self-consistent but WRONG answer here, because the code itself is
wrong — I concluded "Angstrom" from the symmetric import/export and had to be corrected by the A/B.
Do not trust the round-trip as evidence of the unit.

**How to apply:** divide Angstrom shifts by smpd before writing the star. Two other import traps
from the same session: SIMPLE needs RELION **3.1** (`data_optics` + `data_particles`, else zero
orientations import), and put `_rlnImagePixelSize` on EVERY PARTICLE — the optics-group join left
`ogid=0` on the stk records, so `smpd` was missing and `get_ctfparams` hard-failed.
See [[remote-build-recipe-fwl-c143766]].
