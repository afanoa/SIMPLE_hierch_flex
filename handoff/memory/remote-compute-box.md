---
name: remote-compute-box
description: "The reachable remote box (fwl-c142339, 80 cores) — how to connect, what data is there, and what is NOT"
metadata: 
  node_type: memory
  type: reference
  originSessionId: e31f34db-3e71-4d09-ab25-71f10f48fa4c
  modified: 2026-08-07T14:32:37.823Z
---

`ssh atifao@fwl-c142339.ncifcrf.gov` — key auth already works, no port to open.
**80 cores.** Repo at `/home/atifao/SIMPLE`, build dir present and working
(`build/production/simple_exec`); `make -j40` is fine.

Data lives under `/usr/local/data/atifao/`:
- `Ribosembly/` — raw benchmark: `images/`, `ribo_fwd.star`, `ribo_rev.star`,
  **`gt_latents.pkl`** (the ground truth), `vols/`, `pdbs/`
- `IgG-RL/` — the IgG heterogeneity set, and **its project IS built** (unlike Ribosembly):
  - projfile `IgG-RL/igg/4_reconstruct3D/igg.simple`, vol1 `.../4_reconstruct3D/recvol_state01.mrc`
  - 100k particles, box 128, smpd 3.0; stacks are ABSOLUTE paths under `IgG-RL/images/snr0.01/`
  - GT volumes `IgG-RL/vols/128_org/*.mrc` (100 of them); `analysis/gt_frame_coords.npz` is
    regenerable from them by `gt_coords.py <volsdir> <out.npz>`
  - a good reference run to score against: `igg_aug5cfg/FIXED/` (d=128, caps ON, 5:09 wall)
- Ready-made run scripts `igg_*.sh` show the working env: GCC module
  `/mnt/nasapps/development/modules_elmlund/gcc/16.1.0`, `SIMPLE_QSYS=local`, `mkdir=no` with a
  dedicated output dir. A fresh worktree builds clean with cmake Release + USE_OPENMP + gcc 16.1.

**No `mrcfile` in any python on this box** (system, anaconda3). Don't install it — MRC is a
1024-byte header plus NSYMBT then raw data, so a 10-line numpy reader works; `score_states.py`
and `gt_coords.py` both need that swap, and their GT paths are hardcoded to `/home/afan/...`.
- `ribo_setup.sh` / `ribo_run.sh` / `make_ribo_star.py` — project construction scripts
- `tomotwin_flex/env.sh` — environment sourced by those scripts

**`ribo_flex/` is EMPTY** — the Ribosembly *project* is not built on this box.
`ribo_setup.sh` creates it: new_project -> import_particles from `ribo_fwd.star`
-> reconstruct3D consensus (that consensus is `vol1` for flex_pca). Budget an hour
before flex_pca can start.

**This is NOT the machine the 2026-08-06 speedup handoff came from.** That one used
`/home/afan/SIMPLE` and `/home/afan/simple_ribosembly` on a 32-core box, and no
`/home/afan` exists here. Numbers from that handoff are not reproducible here without
rebuilding the project first.

See [[igg-test-project-setup]] for the IgG side, [[flex-pca-speedup-verified]] for what
was validated.
