---
name: rnasep-bins-overnight-2026-10-15
description: where the 2026-10-15 overnight RNase P direction-bin chain lives (batch2 rnasep_coupled/bins80 + clusters; Mac scratchpad overnight3.log) and what it produces (sync windows bsyn1-3, grid bsg, pooled grid bgrid9, montage); the box-160 bins/ results and the first bins80 pass are VOID (scrambled import); read overnight3.log first after a compaction
metadata:
  type: project
---

Overnight 2026-10-15 (user: "work autonomously overnight and see if you can get anywhere"): the direction-bin route
on RNase P. batch2: `rnasep_coupled/bins` (project with class = bin, 100 bins, refined particles only via selection
with flags_refined.txt on ptcl2D+ptcl3D; make_cavgs; cls_expansion pose mean 16 sub-classes per bin; dumps job
`run_bins_dumps.sbatch`; pooled-latent grid `bins/grid/maps` = bgrid9), `rnasep_coupled/clusters` (keep files and
reconstructions `run_grid_rec.sbatch` with ARMS/NG; maps bsyn1-3 (8 windows per synchronised axis), bsg (3 x 3 grid),
grid9 (pooled 2D-class grid, diffuse)). Mac scratchpad: `bins/` (dumps, latent dumps, sync/, win/keep, maps/),
`overnight.log` (the chain's log; look for OVERNIGHT SYNC DONE), figures copied to doc/implementation_notes/
rnasep_coupled_figs/viz_rnasep_bins_*.png. Scripts in scripts/cls_expansion_proto/subclass_pool/ (overnight_sync.sh,
sync_windows.py, windows_readout.py, grid_readout.py, grid_keep.py, bins_montage.py, run_rnasep_bins.sbatch).

**Why:** the flex-latent grid on RNase P gave reproducible but diffuse differences (section 21); the 2D side on
direction bins, stitched by synchronisation, is the remaining route to visible domain movements.

**How to apply (updated 04:40):** the live chain is overnight3.sh (log scratchpad/overnight3.log, markers OVERNIGHT3
DONE/FAILED) on rnasep_coupled/bins80 (row-ordered bin-frame stacks, run_bins_transform.sbatch v2 ->
run_rnasep_bins80.sbatch -> run_bins80_dumps.sbatch -> sync -> run_grid_rec.sbatch ARMS='bsyn1 bsyn2 bsyn3'/bsg);
local results in scratchpad/bins80, figures viz_rnasep_bins80_*.png. Everything produced from rnasep_coupled/bins
(box 160) and from the first bins80 pass (before 04:25) is void: see
[[rnasep-692k-project-rows-are-a-permutation-of-stack-order]]. The bins80 job is a 13-minute loop (import 4 min,
cls_expansion 4 min, nine reconstructions 4 min). Related: [[toyl-direction-bin-scheme-results-2026-10-14]],
[[rnasep-sync-od-coordinate-2026-10-14]].
