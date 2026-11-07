---
name: empiar-10330-cs-corrected-optics
description: "EMPIAR-10330 / PfCRT was collected on a Cs-CORRECTED Titan Krios: nominal_cs = 0.001 mm, not 2.7. SIMPLE's cs default is 2.7 and nothing validates it, so any import that omits cs=0.001 puts ~4 rad of spurious phase into the CTF at 3 A. Plus the deposit's own layout traps (Patch2/Patch3 duplicate basenames, EMDB defocus field is 1000x off, counting mode not super-res)."
metadata:
  type: reference
---

**Authoritative source:** EMD-20806 deposition record (`https://www.ebi.ac.uk/emdb/api/entry/EMD-20806`,
`structure_determination_list -> microscopy_list`), cross-referenced from EMPIAR-10330. The paper
(10.1038/s41586-019-1795-x, open access at PMC6911266) says "Titan Krios electron microscope (FEI) equipped with a
Cs corrector" but never states a Cs value -- only the EMDB record does.

Verified parameters:
- microscope FEI TITAN KRIOS, `sph_aberration_corrector` = "Cs corrector was used.", **`nominal_cs` = 0.001 mm**
- 300 kV; GIF BioQuantum energy filter, 20 eV slit
- Gatan K2 Quantum, detector_mode **COUNTING** (NOT super-resolution): 3838 x 3710 at 0.5175 A/pix physical,
  calibrated on an AAV2 dataset. Nyquist of the raw micrographs is 1.035 A. Cross-check: 4.0 e/pix/s / 0.5175^2
  x 6.0 s = 89.6 e/A^2, matching the stated 91.56 -- the dose arithmetic only closes with 0.5175 as the PHYSICAL pixel.
- 80 frames, 6.0 s, 91.56 e/A^2; defocus -1.2 to -1.8 um; C1; 3.3 A (3.2 A for PfCRT alone)
- PICKING CASCADE (from the paper's methods -- do NOT read 183,241 as the initial pick count, that error drove a
  wrong 'we are over-picking' conclusion once): DoG picker 682,520 (207/mic) -> RELION 2.1 2D -> good classes used as
  Gautomatch templates with a DELIBERATELY LENIENT threshold giving **1,095,304 initial picks (332/mic)** -> pared to
  183,241 after 2D (16.7% kept) -> 17,030 used in the final reconstruction (1.5% of initial picks). Particles were
  extracted BINNED BY 2 (0.5175 -> 1.035), matching our binning. Final map reconstructed in cisTEM from non-subtracted
  particles; angles from cryoSPARC non-uniform refinement.
  Practical consequence: lenient picking + aggressive 2D rejection is the proven recipe here. A SIMPLE refpick at the
  default separation gives ~223/mic, already 33% BELOW their lenient pick, so do not tighten it.
- amplitude contrast never stated in the paper; both deposited star files carry 0.1.

**Why it matters:** the Cs term of the CTF phase is (pi/2) Cs lambda^3 k^4. At 300 kV (lambda = 0.0197 A), using
cs=2.7 instead of 0.001 adds ~0.03 rad at 10 A, ~0.25 rad at 6 A, ~1.3 rad at 4 A and ~4.0 rad at 3 A. So a wrong Cs
is nearly harmless for heterogeneity work in the 8-30 A band but fatal for anything refining past ~6 A. SIMPLE's
default is `cs=2.7` (src/main/params/simple_parameters.f90) and there is no lower-bound validation, so the error is
silent. Worth checking against the 3.93 A refine3D in [[pfcrt-poses-are-the-blocker]].

**How to apply:** import with `smpd=1.035 kv=300 cs=0.001 fraca=0.1 ctf=yes` (smpd 0.5175 if unbinned). Deposit traps:
- `02_Aligned_Micrographs/` holds `Patch3_micrographs/` (3297, 3x3 patches, complete quality-filtered set) and
  `Patch2_micrographs/` (3174, 2x2 patches, a re-alignment of a SUBSET). Same basenames in both -- globbing
  `Patch*_micrographs/*.mrc` imports ~3174 duplicates and `add_movies` will not catch it (it only checks dimensions).
  The deposited `Patch3_micrographs_ctf.star` matches Patch3; the final particle star references Patch2 paths.
- EMDB `nominal_defocus_min/max` read 0.0012/0.0018 um -- a depositor unit slip (mm typed into a um field). The real
  range is 1.2-1.8 um, confirmed by the CTFFIND values in the deposited star (1.07-1.88 um).
- EMDB `nominal_magnification` 215,000x does not reconcile with 0.5175 A/pix (that implies ~96,600x). Irrelevant for
  SIMPLE, which takes smpd directly -- do not derive pixel size from the magnification.
- Micrographs are dose-weighted (`-DW.mrc`) and magnification-anisotropy corrected; no non-DW version is deposited,
  and the depositors fit CTF on the DW micrographs too.
- SIMPLE has NO magnification-anisotropy correction anywhere in src/, and this data needs a 1.3% correction
  (angle 172.7 deg, scales 1.007/0.993) -- so reprocessing from the raw frames inside SIMPLE alone would be worse
  than starting from the deposited aligned micrographs.

**Built project (2026-09-13, batch2.ncifcrf.gov):** `/mnt/beegfs/elmlund/afan/pfcrt_from_scratch/`
- download lives at `/mnt/beegfs/elmlund/EMPIAR_downloads/10330` (lowercase `downloads`). It arrived with ONE truncated
  Patch3 file (4096 B, `..._00009gr_00047sq_v02_00002hln_00002enn-b-DW.mrc`), re-fetched from the EBI https mirror.
  Verify a download by file-size histogram, NOT by `du`: BeeGFS under-reports allocated blocks for freshly written
  files (reported 122M for 47 GB of real data), so `du` vs apparent size looks like sparse holes when nothing is wrong.
  The real test is reading bytes -- a genuine micrograph has ~1 exact zero in 14.2M samples.
- binned 2x with `relion_image_handler --angpix 0.5175 --rescale_angpix 1.035` (module relion/4.0.1; SIMPLE's own
  `prg=scale filetab=` must NOT be used, its writer opens square and these micrographs are rectangular). 3838x3710 ->
  **1918x1854** (3838/2 is odd so RELION rounds to even), so the true pixel size is **smpd=1.03554**, not 1.035.
  10-task slurm array x 8 procs binned all 3297 in 2 min.
- **Particle stacks inherit the micrograph smpd exactly** -- `extract` is pure windowing, opening the output stack with
  `params%smpd` (simple_extract_strategy.f90). So stacks will also be 1.03554, vs the deposited stack's 1.035.
- imported with `simple_exec prg=import_movies filetab=... smpd=1.03554 kv=300 cs=0.001 fraca=0.1 ctf=yes`; project at
  `pfcrt/1_import_movies/pfcrt.simple`, 3297 records, `imgkind=mic`. Binary `/home/atifao/SIMPLE/build` (f79d1402).
- `script=yes` and `qsys_*` are INERT for import_movies: `script` is only read inside qsys_env, and import_movies
  never builds one. They do apply to ctf_estimate. The key is `qsys_name`, not `qsys`. The project's compenv segment
  already came out as `qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=16000`.

**ctf_estimate (2026-09-13):** `nparts=10 nthr=8 dfmin=0.8 dfmax=2.5 astigtol=0.3 lp=3.5` (hp=30, pspecsz=512,
ctfpatch=yes defaults). 3297 micrographs in 736 s over 10 parts. VALIDATED against the deposited
`Patch3_micrographs_ctf.star`: median |mean-defocus difference| = 0.0022 um (22 A), 97.6% within 0.05 um,
r = 0.945 overall and **0.995 excluding 25 outlier micrographs (0.76%)**. This is the gate that proves cs=0.001,
smpd=1.03554 and the binning are all correct -- run it on any re-import.
- `dfmin=0.8` is safe: only ONE micrograph in the whole set sits below 0.85 um.
- Convention note: SIMPLE dfx/dfy are swapped relative to RELION DefocusU/V and angast differs by ~90 deg. Compare
  MEAN defocus, not U-to-dfx.
- SIMPLE astigmatism runs lower than CTFFIND's (median 0.069 vs 0.108 um) though correlated at r=0.86.
- ctfres median 3.28 A, floored by our lp=3.5 (binned Nyquist is 2.07); deposited got 2.54 A on UNBINNED data.
  192 micrographs fit worse than 5 A -- rejection candidates before picking.

**cluster_cavgs / select_clusters (2026-09-13, general SIMPLE workflow):**
- `cluster_cavgs` (affinity propagation, 43 s on 200 classes) writes per-class and per-particle `cluster` labels into
  the project, ranks clusters by `jointscore`, and **also writes its verdict**: particles in clusters it judges bad get
  `state=0`. `rank_clusters` renumbers `labels` so cluster index 1 = best, so file `clusterNN_cavgs.mrc` DOES match the
  log's `cluster_rankedNN` line (the log names a file that is never written -- don't go looking for it).
- **To override its verdict, select on the `cluster` label ALONE.** Filtering on `cluster in {...} AND state==1` returns
  only the clusters SIMPLE already liked -- on 10330 that gave 23,049 instead of the correct 162,168.
- **Use `prg=select_clusters projfile=... clustinds=4,6,7,8,...`** -- a one-liner that does the whole thing. It PRUNES
  (the project keeps only the selected particles, not all entries with state flags). Verified identical to a hand-built
  `prg=selection oritype=ptcl2D infile=<0/1 per entry>` (both 162,168). Don't hand-roll the 0/1 file.
- **jointscore is a bad "is this a real particle" metric on small/low-contrast specimens.** It is built from resolution
  + homogeneity + compactness, and resolution rewards compact high-contrast averages -- so aggregates and micelle-only
  blobs score ABOVE genuine views where a Fab protrudes into noisy solvent. On 10330 the auto-KEEP set (clusters 1-5)
  was mostly aggregate/micelle and the real particles sat in the REJECTED clusters. Same reason per-class `rank_cavgs`
  ordering was uninformative there (ranks 3-168 all 8-11 A). Always eyeball the cluster montage.

