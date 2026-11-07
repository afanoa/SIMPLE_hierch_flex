---
name: pre-release4-projects-fail-flex-contracts-2026-10-22
description: "TRAP: SIMPLE projects written before release 4 (2026-10-03; e.g. cnga1 flex_pcg_env/proj.simple, PfCRT) stop current flex_pca: ptcl2D/ptcl3D eo labels differ, and the stack lacks nptcls_stk; fix_projfile refuses when the stack holds more images than rows (no indstk contract) -> re-import or re-extract"
metadata:
  type: project
---

2026-10-22, cnga1 complement test. The September cnga1 project (245,424 particles from final_particles.mrcs, 274,333 images) failed in
two steps on the current code:
1. "inconsistent even/odd partition" (eo_transaction): ptcl2D and ptcl3D carry independent eo splits. A copy with ptcl2D eo = ptcl3D eo
   passes (binary particle records, slot I_EO = 11).
2. "nptcls_stk not present in stack field" (map_ptcl_ind2stk_ind). `prg=fix_projfile` fills nptcls_stk from the stack header but then
   REFUSES: rows < images and "no recorded indstk contract; the physical image of each particle cannot be recovered (re-extract or
   re-import these particles)". The old indstk values (3..274332, increasing) are not trusted.

A hand byte-patch of the stack record was denied by the permission check. Do not retry it.

**How to apply:**
- Before planning a run on an old dataset, check its project's date against release 4.
- Run `prg=fix_projfile` on a copy first.
- If it refuses, the dataset needs a re-import (or re-extraction) before any current-code run. Budget that time up front.
- No cnga1 or PfCRT project newer than 2026-10-03 existed on batch2 at this date.

Related: [[simple-star-pose-import-traps]], [[import-particles-resets-state-and-eo]].
