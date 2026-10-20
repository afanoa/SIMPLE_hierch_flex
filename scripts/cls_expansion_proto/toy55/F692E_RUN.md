# f692e snapshot: RNase P real data, delivered 2026-10-18 14:05

This commit preserves the code that produced the f692e results on the 692k refined RNase P particles:
63 leaf states, 212 zoomed sub-states including the pure closed state (zoom state 1, closed-arm occupancy 0.94,
4,896 particles), and the composed maps on per-family parents.

## Contents

- **klem board (Fortran):** `src/main/commanders/simple/simple_klem_{core,distr,tree}.f90` and
  `simple_commanders_klem{,_board,_tree}.f90`, plus the three SIMPLE files the klem2 worktree modified
  (`production/simple_private_exec_driver.f90`, `src/main/exec/simple_exec_refine3D.f90`,
  `src/main/ui/simple/simple_ui_refine3D.f90`) and `register_klem.py`. Base commit 5e1e6aeee, the same base as
  the batch2 worktree `~/SIMPLE_klem2`. These sources match the binaries the run used, which were copied on
  2026-10-17 18:36 to `/mnt/beegfs/elmlund/afan/rnasep_klem/f692e/klem2_bin` on batch2.
- **Pipeline scripts:** `scripts/cls_expansion_proto/toy55/` as of 2026-10-18 14:05. The klem sources are mirrored
  in `toy55/port/`.

## How the run was made (batch2, `/mnt/beegfs/elmlund/afan/routes/toy55`)

```
sbatch --export=ALL,RTAG=f692e,DEV=0,LOCALDEV=0,SERVER=1,NPARTS=20,PARMAX=4,CTFABS=1,PSCALE=1,REFINE_EXTRA=ml_reg=no,NFREE_SUB=1,NP=20 f692_chain9.sbatch
```

Then, on the finished run: `zoom_leaves.sbatch` (RTAG=f692e, FAM=1 and FAM=2), `zoom_analysis.sbatch`,
`pooled_dd.sbatch` (LEAFGLOB="region*" FLOOR=12 TAG=dd, and LEAFGLOB="zoomregion*" TAG=ddzoom), the all-particle
consolidation (`consolidate.sbatch` with UNCAP=1, giving `consolidated_12_fullw_region` and
`consolidated_12_fullw_zoom`), and `compose_state_maps.py` with per-family parents at the families' refined poses.

Outputs: `/mnt/beegfs/elmlund/afan/rnasep_klem/f692e`. Delivered copies:
`~/Downloads/rnasep_klem_20261015/deposited_conformers/{f692e_leaf_states_allptcls,f692e_closed_substates,f692e_composed_maps}`.

## How this snapshot was reconstructed

None of this code was under version control, so later edits were undone from the recorded edit commands.
Every undo was checked by re-applying the original edit, which reproduces the 2026-10-19 file byte for byte.

- `simple_klem_core.f90`, `simple_klem_distr.f90`, `simple_commanders_klem_board.f90`: three later additions
  removed, all off by default: Kaiser-Bessel gridding (2026-10-18 10:52-11:37), `specmatch` (2026-10-19 10:23)
  and `hp` (2026-10-19 12:49).
- Five scripts changed after the cutoff were restored: `f692_family.sh` (the `REFINE_EXTRA` default later became
  `filt_mode=none`), `consolidate.sbatch`, `consolidate_states.py`, `merge_data_derived.py` (the merge rule
  changed after this run) and `pooled_dd.sbatch`. The f692e family refinement log shows the default nonuniform
  filtering, which confirms that `filt_mode=none` was added later.
- 167 scripts created or changed after the cutoff are not in this commit. The next commit on this branch adds the
  state of 2026-10-19.
