---
name: rnasep-full81k-preprocessing-run-2026-10-22
description: "Full RNase P apo reprocessing from movies (2026-10-22): 5 grids, 81,265 EER movies read in place, two gain groups (no flip), smpd 0.732 -> 1.464, dose 53, pick refs = top-40 ranked 4_abinitio2D cavgs (user: includes the closed state), box 160; /mnt/beegfs/elmlund/afan/RNAseP_Full"
metadata:
  node_type: memory
  type: project
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-10-24T03:00:14.445Z
---

Working dir `/mnt/beegfs/elmlund/afan/RNAseP_Full` (raw data only read, never copied). Build `~/SIMPLE_proc` (aea5b1846).
- Grids (all Krios-G4 9930309): Nov 105640_5_apo 11,360, grid5 18,777, grid3 10,056, gri10 9,644 (gain 20241031) and
  Dec 4_2_b6 31,428 (gain 20241210, a DIFFERENT file) -> projects `nov` and `dec` (gainref is per run), merge later.
- Gain orientation tested on both (fixed-pattern ratio R): flipgain=no is right for both.
- EPU XML pixel is 0.7048 (nominal); use the calibrated 0.732 (paper). Dose 53 e/A2 from the paper's SI Table 1
  (grid5/Dec XMLs lack a dose field). preprocess smpd=0.732 smpd_downscale=1.464: frames Fourier-cropped BEFORE alignment,
  CTF fit on the 1.464 forctf sum (30-5 A).
- Pick refs: stackops top 40 of `afan/atifao/reextract/4_abinitio2D/cavgs_iter044_ranked.mrc` (box 160 @ 1.464);
  the user checked the montage: the closed state is among them.
- Pilot (250 movies, `pilot/`): 11.5-11.9 s/movie per 48-core node but only ~1 core busy per part (per-part bound) ->
  run many 4-thread parts; ~20 GB per part; 271-282 picks/mic, ~92 % extracted; ctfres 3.7-4.7 A, gri10 6.2 A (worst).
- Full run: run_full.sbatch (import -> preprocess -> pick -> extract box 160), workers on norm at the coordinator's
  request (dec 47 parts 54106937, nov 75 parts 54108982), 4 thr / 24 GB / time_per_image 150 (workers got 27.7 h).
  SIMPLE sets worker --time = time_per_image x partsz capped at 48 h; norm/csbdevel allow 5 days. User 2026-10-23:
  raise it -> `export SIMPLE_DEFAULT_PARTITION_TIME=true` (no --time line, partition default 5 d); run_full.sbatch
  now does this (old copy run_full_v1.sbatch).
- DONE 2026-10-23 18:54: preprocess ~17.4 h per group (slow tail), pick 80-100 min, box-240 extract 14-17 min.
  19,280,495 particles at box 240 / 1.464 A from 81,265 mics (237/mic): dec/5_extract 7,347,311 and nov/5_extract
  11,933,184 (projects dec.simple / nov.simple there, not yet merged). dec|nov/4_extract = cancelled box-160 stub only.
- Extraction switched to box 240 (351 A, the paper's classification box): guard jobs cancel the masters' box-160 extract
  after pick; run_extract240.sbatch extracts the picks again at 240. Masks 240-250 A in 2D and 3D.
- USER PLAN (2026-10-23, via coordinator, revised the same day): expect ~4-5 M clean particles; the full-data tree runs on
  ALL clean particles (~12 h wall accepted, both partitions), NOT subset-then-assign. Still keep every 2D reject (stacks +
  metadata) findable for possible later use: rejection = state 0 only, never prune_project in place; write explicit row
  lists clean / rejected / all and report the counts to the coordinator.
- 2D cleaning 2026-10-23: full/ = merge (1_merge_projects, rnasep_full.simple) -> rnasep_full_q.simple (compenv norm,
  16 GB, 4 h) -> particle_sieving single_pass=yes mskdiam=240 50x8 (job 54180123, 3,767 chunks, ~79 % kept). Row lists
  -> lists/. Waiter 54182549 then submits run_solve2d.sbatch (solve2D 500 cls nsample 400k mskdiam 240 48x8 + cluster_cavgs).
  User: "its fine you dont have to ask" = launch solve2D without the price check (still report the price).
  update_project projfile_out needs a root mapping: copy the project, then update in place.
- SIMPLE's slurm scripts hard-code `--mail-type=FAIL` (simple_qsys_slurm.f90); mail-user is nil.
Related: [[rnasep-degenhardt-paper-closed-state-and-pipeline]], [[never-delete-only-read-or-move]],
[[rnasep-closed-state-not-pure-in-klem-leaves]].
- Time limits (user 2026-10-25, "why is the time limit 18h thats stupid"): no tight -t from "2x the estimate" rules; masters
  -t 5-00:00:00 and workers SIMPLE_DEFAULT_PARTITION_TIME=true (partition default, 5 days on norm/csbdevel).
- Cost accounting: solve2D workers hold cores only ~3 min per ~16 min iteration (master-bound on the 10 GB project), so
  price from sacct worker core-h, not cores x wall: round 1 = 1,971 worker core-h + ~120 master, not 5,900.
- Full-res test 2026-10-25 (full/9_class1_fullres, user asked): round-2 rank-1 class 70, the same 2,990 particles
  re-motion-corrected at 0.732 (box 480) vs 1.464 (box 240): binned 3.70 A (SIMPLE) / half-FRC 0.5 at 5.2 A; full 4.7-4.8 A /
  6.8 A, no signal beyond 3 A; same at matched 3.0 A alignment limit. No gain from full res in 2D at this N; keep 1.464 for
  classification (as the paper did); revisit full res only for final per-conformer 3D. Cost 488 core-h, ~1.5 h.
