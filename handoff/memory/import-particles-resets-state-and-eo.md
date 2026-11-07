---
name: import-particles-resets-state-and-eo
description: TRAP: import_particles (oritab + stktab) sets state=1 on EVERY row and eo=0 regardless of the oritab; a reconstruct3D on the imported project then uses all particles (the "10k" re-imports were 932k maps, FSC 4.3 A, arm = family average); run prg=selection with a stack-order states file after the import
metadata:
  type: project
---

2026-10-16 02:45: three re-imports of the RNase P within-10k set (ctf=yes, shifts flipped, ctf=flip) all gave contact
0.30, arm 0.42; pairing, poses and stacks were verified identical to the family project. The imported project's
ptcl3D field had ` state=1` on all 932,496 rows and ` eo=0` everywhere: the oritab's state column is ignored.

**How to apply:** after `import_particles ... oritab=... stktab=...`, always
`simple_exec prg=selection projfile=... oritype=ptcl3D infile=<states in STACK order>` and verify with
`print_project_field | grep -c ' state=1\.'` (note the leading space: `mi_state=` also contains "state=").
Pass `eo=` in the oritab or check how reconstruct3D splits halves. The stacks are still `ctf=flip`
([[rnasep-denoised-stacks-are-phase-flipped]]). The fixed job is region_align_rec_rnasep.sbatch (unmoved control +
region + core2).

- 2026-11-06 refinement (fixture session): the reset was seen with `stktab=` imports. A single `stk=` import plus an oritab
  with `eo=` kept eo balanced (IgM fixture 50,000 / 50,000). Still print the state=1 / eo=0 / eo=1 counts after every import.
