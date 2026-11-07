---
name: klem-f692-fortran-chain
description: the unsupervised klem tree on the full 692k apo set with the Fortran board runs as ONE resumable sbatch chain (toy55/f692_chain.sbatch -> fboard.sh, f692_family.sh, f692_regions.sh); outputs under rnasep_klem/f692; traps (refine3D_states maxits >= 5, frozen build copy, smoke mode)
metadata:
  type: project
---

2026-10-16: `f692_chain.sbatch` (toy55, mirrored on batch2) = parent Fortran board on all 692k at the consensus poses
(7x7 + 2 free, 20 -> 10 A, latent init from meta_all z692) -> expand_families.py split (2 families) + derive_regions
on the parent -> per family `f692_family.sh` (child board, k-means K=10, refine3D_states local 20 -> 8 A maxits 6
nparts 10, feedback board at the refined poses with the child's own region as focus) -> per family `f692_regions.sh`
(derive_regions + derive_regions2 on the feedback board, one Fortran focal board per region, node_states, recs,
level 2 ranked along the node's own axis) -> count_states over the leaves. Every stage resumes from existing outputs
(resubmit with the same RTAG). `fboard.sh OUT PROJ SUBSET "LPS" "SIGMAS"` wraps prg=klem_board with the validated
conventions (conv 0 ss -1 fs 1 sign 1 ctfabs 0, shrink 0.02, eps0 0.05, balance 0.5, devw 0.05, pnoise 1); env FB
(build dir), META (latent init), N1 N2 NFREE BAL FOCUSMASK INIT POSE_EVERY/STEP/PRIOR. The chain copies the klem2
build to $R/klem2_bin so the worktree can be rebuilt while it runs. Submitted: full run 53905193 (RTAG f692), smoke
(SMOKE=1: 50k subset, short schedules) 53905192 / 53905211. TRAP: refine3D_states refuses maxits < 5 (the first smoke
failed on maxits 1; f692_chain2.sbatch has MAXITS=5). Fortran board memory: 409k at 8 A = 8.5 GB MaxRSS (Python 112 GB).
Not in the chain yet: the docked-node state enumeration + pose polish (f692_docked.sh, run after the chain).
Related: [[klem-fortran-port-state-2026-10-16]], [[rnasep-klem-real-data-setup]].

NEGATIVE (2026-10-16 12:30, toy "ah" with truth, 200k particles, 18 conformers + junk): the Fortran board's whole-reference
local pose step at FAMILY scale (7x7 mixed-state board, pose_every 2, step 3, prior 4) does not replace refine3D_states:
from 1.4 deg poses it worsens them to 3.8 deg, from a realistic 3.8 deg start it gives 3.9 deg (p90 5.5 -> 6.9). It pays
only on pure sets (wide toy 3.8 -> 2.9, real docked 10k). The chain keeps refine3D_states (~3.5 h over two families).
Also: the tied two-layer board (--tie_layers --tie_support) collapsed on the toy (component -> a 0.2 % spike, free comps
took 99 % of the particles): the component needs amplitude clipping to the base density and a harder shrinkage; deferred
(the plain lattice already separates holo from apo on the deposited-holo toy: 21 cells, 74 % of holo at 82 % purity).
Run-time controls added: $R/maxreg_fam{1,2}.txt override the regions per family (set to 4 and 2 to finish ~21:00).

13:05 PARALLEL REGIONS: f692_regions.sh on batch2 is now the parallel version (f692_regions_par.sh: one
f692_region_one.sbatch job per region, PARMAX concurrent, DONE markers, resumable; the sequential original is
f692_regions_seq.sh). Run-time override files in $R: maxreg_fam{1,2}.txt, parmax.txt, nthr.txt (current run: 4/2
regions, 2 concurrent, 16 threads -> finish ~19:30). Next-run chain = f692_chain3.sbatch (-c 32, env NTHR 32, NP
refine3D_states parts (20 at night), PARMAX 3) ~5.5 h; refine3D_states (3.3 h of it) is the floor, the board is not.
The distributed board (workers exchanging partial Fourier sums per round) was costed at a day for ~3x per board and
deferred: with parallel regions and 32 threads it would save only ~30 min of the chain.

TWO TRAPS 2026-10-16 afternoon: (1) distributed-board workers (and tree sub-jobs) must run through simple_private_exec:
simple_exec appends a job record to the PROJECT FILE at exit (update_job_descriptions_in_project), and several workers
exiting together corrupted the toy project's compenv segment (qsys_name lost -> "class: unsupported in qsys_factory").
Repair = prg=update_project with the qsys keys again. klem_board_part and klem_tree are now registered in the private
driver (register_klem.py) and launched with the default executable. (2) NAME CLASH: the regions stage exported N2 (the
level-2 particle count, 10000) and fboard.sh reads N2 as the lattice's second dimension -> 7 x 10000 boards, 115 GB,
OOM. Fixed: the region job renames it NTOP2 and unsets N2. Backups of the real refined projects: *.bak_1430.
Fortran tree driver: port/simple_klem_tree.f90 (families, regions, states, pc1 ranking, identification, landscape) +
port/simple_commanders_klem_tree.f90 (prg=klem_tree, stages all/family/region/node as sub-jobs, config file); builds;
toy validation = run_ftree_toy.sbatch (RTAG toyah_ftree) vs the bash chain toyah_tree_nr.
