---
name: klem-board-rewrites-shared-project-trap
description: "TRAP: every simple_exec prg=klem_board run rewrites its projfile (appends a jobproc record); concurrent boards on one shared project corrupted its compenv -> later refine3D_states on copies die 'class: unsupported in qsys_factory' -> node_refine falls back SILENTLY; run boards on private copies (dtree_board7)"
metadata:
  type: project
---

Found 2026-10-20 02:00: synth_55_700k_B/src.simple had ~194 klem_board jobproc records and its compenv segment held a
klem_board job record (prg=klem_board projfile=... infile=.../sixc8c/r.4/board_s1/board.txt) instead of qsys_name etc.
(concurrent writers interleaved segments). ptcl3D was intact (md5 of e1/e2/e3/x/y/dfx = the pose table), so board and
bench results stand; but every refine3D_states on a copy failed in 3 s ("ERROR! class: unsupported in qsys_factory
constructor") and node_refine/oracle_refine FELL BACK to the input poses: sixc8c r.4 B, oracle B, polish B, all of
sixc8d/sixc8s (v14; those also hit filt_mode=uniform), real gs2 r.3.4 A (depth-2 boards rewrite the shared depth-1
refined project). Check: grep FAILED in node_*.txt.ref; compenv via print_project_field oritype=compenv (compenv_check.sbatch).
Real data: synth_rnasep692_A/B src.simple clean; gs2 r.3 A refined project corrupted (72 board records) with ptcl3D intact
(projcheck.sbatch md5) -> only r.3.4 A lost its depth-2 refinement. Repair: compenv_repair.sbatch (update_project qsys_name=slurm ...; keeps a .bak_compenv; stray keys remain).
Fix: dtree_board7.sbatch copies PROJ to the board dir and deletes it after. Same risk for any program run concurrently
on one projfile (node_refine already copies). refine3D_states takes filt_mode=fsc|nonuniform_lpset|none only.
Related: [[klem-worker-scripts-collide-in-shared-cwd]], [[klem-data-derived-schedule-2026-10-20]].
