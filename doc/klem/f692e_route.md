# The f692e klem route as one SIMPLE program

`prg=klem_tree` runs the kernel-latent EM route that produced the f692e RNase P landscape (63 leaf states, the pure
closed sub-state, the composed maps) and its rerun f692r. Before this, the route was ~14 bash/sbatch scripts and 18
Python helpers in `scripts/cls_expansion_proto/toy55/`; tag `f692e-2026-10-18` keeps that code exactly as it ran.

## Running it

```
sbatch --export=ALL,KLEM_BUILD=$HOME/SIMPLE_klemc/build,PROJ=/path/rnasep.simple,CFG=scripts/klem/rnasep_f692.txt scripts/klem/klem_tree.sbatch
```

`CFG` holds `key value` lines (`klem_tree_cfg` in `simple_commanders_klem_tree.f90`); the defaults are the f692e
settings. `run` is the output directory, `meta` the npz with the initial latent (`z692`) and the particle flags,
`cons`/`docked`/`docked_all`/`closed` the consensus, the docked deposited maps and the closed conformer (needed for the
merge, readouts and the arm readout), `subset` an optional parent subset (one 1-based project row per line), `partition`
an optional SLURM partition for every job. Every stage resumes from its outputs: resubmit the same command after a
failure. The launcher refuses a non-empty `run` directory that it did not start, so a run never resumes from another
route's outputs. A run keeps its private copy of the binaries in `run/klem_bin`; remove `run/klem_bin/bin` before
resubmitting to resume on a newer build.

## Stages

Each stage is the same program with `stage` set; the chain submits the others as SLURM jobs through
`qsys_async_job`, so a job that crashes or is killed fails the run instead of being waited for.

| stage | what it does | bash/Python it replaces |
|---|---|---|
| all | parent subset, parent board, families, parent region, family jobs, landscape count, zoom jobs, pool and parent jobs, composed maps | f692_chain9.sbatch, expand_families.py, derive_regions.py, count_states.py, compose_state_maps.py |
| family | child board, k-means K=10, refine3D_states, feedback board as the family's board server, regions, region jobs | f692_family_job3.sbatch, f692_family.sh, f692_regions_par3.sh, derive_regions2.py |
| region | focal board on the server, supported states, reconstructions, level 2 along the node's own axis | f692_region_one3.sbatch, node_states2.py, rank_by_target.py |
| zoom | one leaf state: zoom board, sub-states, local reconstruction | zoom_leaves.sbatch, zoom_one.sh, zoom_rescue.sbatch |
| pool | data-derived merge of the region or zoom leaves, all-particle consolidation, gallery, arm readout, recovery | pooled_dd.sbatch, merge_data_derived.py, consolidate.sbatch, consolidate_states.py, pooled_analysis.sbatch, volume_gallery.py, arm_readout.py, recover55.py |
| parents | per-family parent maps at the refined poses | parent_fam_compose2.sbatch |

`klem_init_positions.py` (lattice start positions from the npz latent), `klem_fortran2npz.py` (no longer needed),
`deposited_match.py` (identification readouts) and `repro_compare.py` (now `cmp`) are ported as well.

## Behaviour kept from the 2026-10-20 fixes

- Every board that is not served runs on a private copy of its project (`<board>/p.simple`).
- Reconstructions switch off states with fewer than `minrec` (50) particles or fewer than `minhalf` (10) in the even or
  odd half; even/odd is read from the project.
- Zooms run one job per leaf state (`zparmax` at a time); a zoom whose reconstruction failed is redone on resume.
- `subset` replaces SUBSETFILE; `partition` replaces the norm/unl copies of the scripts.
- The map folders of the results (consolidation `maps/`, `composed/`, `parents/`) hold full maps only; their even/odd
  halves are in a `halves/` subfolder.

## Differences from the bash/Python route

- The two k-means steps (family states K=10, motion-pattern groups of the regions) use SIMPLE's random numbers instead of
  numpy's, so they match the Python only statistically. Every other step is a numerical twin of its Python script.
- A failed refine3D_states or reconstruction stops the run; the bash fell back to the consensus poses.
- The gallery is written as JPEG pages (projections and deviation projections) plus the same table.

## Regression test

`scripts/klem/klem_board_repro.sbatch` reruns the Oct 17 smoke parent board (`f692smoke7`) on the frozen f692e
binaries and on `KLEM_BUILD`, and compares every output byte for byte (`VERDICT: IDENTICAL`).
