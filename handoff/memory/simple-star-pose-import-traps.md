---
name: simple-star-pose-import-traps
description: "Importing a RELION star WITH poses into SIMPLE: reader drops each row's LAST column, stk box gets optics' bogus value, and relative image paths resolve against the star's own directory"
metadata:
  type: feedback
---

Validated recipe for getting external per-particle poses (e.g. a cryoSPARC refinement)
into a SIMPLE project, plus the three traps that each cost a debugging round
(EMPIAR-10028 cryoSPARC-pose import, 2026-08-22, probes A-H):

**Recipe:** RELION 3.1 star with `_rlnAngleRot/Tilt/Psi` (deg) + `_rlnOriginX/YAngst`
written in **pixels** (see [[simple-star-shifts-are-pixels]]) + **one trailing pad
column** -> `new_project` -> `import_particles starfile=... dir=... smpd= kv= cs= fraca=`
-> **byte-patch** `box=3.0000000000000000` -> `box=360.00000000000000` (equal length)
in the .simple -> verify ptcl3D e1/e2/e3/x/y against source values. Working converter:
`~/pose_convert.py` on fwl-c142339.

1. **The star reader silently DROPS each data row's final column.** With
   `_rlnOriginYAngst` last, every y imported as 0 while e1/e2/e3/x were exact -- no
   error anywhere. Fix: append any harmless label after it (`_rlnRandomSubset` -- also
   gives an eo split). Check y explicitly after ANY star pose import.
2. **With pose columns present, stk records get `box=3` (the optics segment's bogus
   box) instead of 360 from the stack headers** -- every downstream program dies with
   "box needs to be larger than 26", which reads like a project-wide corruption but is
   1082 text substrings. Same-length byte-patch fixes it; segments store text records.
3. **Relative `_rlnImageName` paths resolve against the STAR FILE's directory**
   (`rootdir = stemname(filename)` in `import_ptcls`), not `dir=`, not cwd. A star
   outside the data tree prunes every particle ("no stack found"), which then causes
   trap 2's box error with a misleading signature. Absolute paths bypass rootdir
   entirely (leading "/" check) -- but only in `import_starproject`'s path;
   `import_particles` wants the star next to `Particles/`.

Also learned: `import_particles` **ignores `oritab=` whenever `starfile=` is given**
(the star branch returns before oritab is read), and `import_starproject
import_type=ptcl3D` fills ONLY the ptcl3D table -- into a fresh project it leaves the
stk segment empty, which is unusable. `import_particles`' own star path imports poses
into ptcl3D correctly (particles3D flag table maps rlnAngleRot->e1 1:1, degrees,
no unit conversion on origins).

**Trap 4 (2026-09-04, cnga1):** a star-imported project LOSES sigma2 discovery
— flex falls back to the unit white-noise spectrum (loud "SIGMA: provide
sigma2_group*.star" line) and on cnga1 that alone collapsed state eo-FSC from
4.4-5.0 A to 17-22 A with a WILD per-shell mean scale (5.4/2.4/2.0/0.04...).
Fix: stage the source refinement's sigma2_it_*.star into the run cwd. The
diagnostic chain that found it: copy0-only round-trip broke identically ->
poses bit-identical -> diff the run logs, the SIGMA fallback line was the
divergence. Also: import_particles DEDUPLICATES star rows by image index
within one stack (981k rows -> 245k silently, keeps last) — symmetry
expansion needs one symlinked stack PER symmetry copy (zero bytes,
final_particles_symK.mrcs) so each copy is a distinct stack; write the star
copy-major and verify with a stride-N check, and C4 expansion itself is just
rot += k*90 (e1 is the innermost molecule-frame z euler: R = Rz(e3)Ry(e2)Rz(e1)).

**Trap 5 (2026-10-06, thyroglobulin):** a BARE `_rlnImageName` (`000001@stack.mrcs`, no `/`)
is never resolved at all -- `simple_starproject` only builds candidate paths when
`substr_ind("/") > 1` -- so it is looked up relative to the project cwd, every particle is
pruned ("no stack found ... pruning referenced particles"), import still exits 0 with an
empty ptcl3D. Fix: put the star one level ABOVE the stack and write `NNN@<stackdir>/stack.mrcs`.
Also: eo=0 on import is fine -- for ptcl oris eo=0 counts as absent, so rec3D/refine3D
call partition_eo on first use.

**Trap 6 (2026-10-06/07, thyroglobulin C2 expansion, master 311b8817a):** SYMLINKED copy stacks no
longer work -- the importer resolves symlinks to the real path (simple_abspath), so every copy becomes
ONE stk record and the dedup keeps one row per image index (67,480 rows in -> 33,740 imported, import
rc=0, NORMAL STOP, 201 s). Use a HARD link per copy instead; BeeGFS only allows hard links within the
SAME directory, so put `stack_symK.mrcs` next to the original stack. Also: flex_pca halves split by
row index mod 4 ({0,1}|{2,3}), so the copy offset N must be == 0 (mod 4) or a particle and its copy land
in opposite half-fits (drop N % 4 particles); the row-alternating eo split then also matches.

