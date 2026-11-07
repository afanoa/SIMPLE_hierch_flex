---
name: rnasep-692k-project-rows-are-a-permutation-of-stack-order
description: "TRAP: SIMPLE reads images via per-row stkind/indstk keys; in rnasep_coupled/full692k_2d/rnasep.simple those keys are a random permutation of the stack order, and import_particles re-keys rows sequentially, so stacks fed to import must be laid out in project ROW order (not by the exported keys)"
metadata:
  node_type: memory
  type: project
  originSessionId: 5b04614c-23e0-4650-8600-03b464a7f576
  modified: 2026-10-15T08:25:09.869Z
---

SIMPLE maps a particle row to its image through the row's own `stkind` and `indstk` keys (`map_ptcl_ind2stk_ind`,
src/main/project/simple_sp_project_ptcl.f90); the os_stk fromp/top ranges are only a fallback. In the RNase P 692k
project ($R/full692k_2d/rnasep.simple, 932,496 rows over 22_denoise_project/denoise_raw_particles_part01-10.mrcs)
the keys are a perfect permutation of the stack order: row 1 is (stack 2, image 30899), row 96140 is (4, 31673).
Row i is NOT image i.

`import_particles stktab= oritab=` ignores the oritab's keys and re-keys rows sequentially (row i -> image i of the
concatenated stacks). So an oritab exported from this project (print_project_field) imported against stacks laid
out by the exported keys pairs every pose with the wrong image: all class averages become one identical blob
(global average, pairwise 0.995), cls_expansion sub-classes become view clusters of random particles ("rotated
copies of the view", a blob bulk sub-class). The first direction-bin passes (box 160 bins/, first bins80/) were
exactly this and are void (2026-10-15).

**Why:** cost a night of the direction-bin scheme; the symptom (identical make_cavgs averages) was first blamed
on the MRC MZ header, which was a real but secondary bug.

**How to apply:** when writing stacks for import from an oritab of this project, place row i's image at sequential
position i (part p = the original fromp..top range) and read the source image via the row's keys
(run_bins_transform.sbatch v2 does this). Reading raw members in Python must use the keys. Any project made by
import from a row-ordered oritab has row i = original row i, so flags_refined.txt, flex coordinates and keep files
over the original rows apply directly. Check make_cavgs bin averages are distinct before trusting anything
downstream. See [[rnasep-bins-need-bin-frame-stacks]], [[rnasep-24-refine3d-auto-quarter-unrefined]].
