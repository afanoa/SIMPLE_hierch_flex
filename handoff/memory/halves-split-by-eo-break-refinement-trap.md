---
name: halves-split-by-eo-break-refinement-trap
description: "TRAP: the real RNase P halves were split by the master project's eo, so a per-half run on the MASTER project has every particle even (A) or odd (B) -> refine3D_states 'canonical sigma2 has an empty even/odd half' in 2 min -> f692 family falls back to consensus poses SILENTLY; run halves on synth_rnasep692_{A,B}/src.simple (import re-draws eo) with a sliced meta"
metadata:
  type: project
---

Hit 2026-10-20 11:08 on the first per-half f692e runs (f692hA/f692hB: SRC = master rnasep.simple + SUBSETFILE half
rows): all four family refinements died ("SIGMA2 GLOBAL BOOTSTRAP ACTIVE E/O: 159698/0", "ERROR! canonical sigma2 has
an empty even/odd half; exec_calc_pspec_assemble") and f692_family.sh logged "no refinement ... at the consensus poses"
and carried on. Zoom reconstructions (reconstruct3D) would hit the same. Fix used: SRC = $K/synth_rnasep692_{A,B}/src.simple
(half projects, mixed eo; per-node and polish refinements ran on them), METAALL = f692e_halfmeta/meta_half_{A,B}.npz
(make_half_meta.py: meta_all per-row arrays sliced by simrnasep692_{A,B}/meta_rows.txt; the chain scripts do not read
meta eo/stkind), no SUBSETFILE; runs f692hA2 / f692hB2 (jobs 54022168/9, norm, LOCALDEV).
Lesson: a fallback that only logs is a silent failure; the job monitor ([[job-monitor-batch2]]) now alerts on
"no refinement|falling back|FAILED".
Second form (2026-10-20 13:46, f692hA3/hB3 on the half projects): family eo balanced (80906/80478) but REGION states
eo-skewed up to 1575/28695 (A fam1 motion_5; B fam1 region_1 7878/22190), never in f692r (all ~50/50) -> a state with
an empty half segfaults reconstruct3D. Guard in routes/toy55_f692e_unl (f692_region_one3.sbatch rec(), zoom_one.sh):
switch off states with < MINHALF (10) particles in either half (eo per row cached in $D/eo_rows.txt via
print_project_field). Cause of the skew not yet explained.
