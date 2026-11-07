---
name: flipqr-dataset-setup
description: "The flipqr test dataset on batch2: path is elmlund NOT emlund, 689,169 particles at box 180 / smpd 1.3 (the 0.822 in the log is the MOVIE smpd), particle ~120 A so mskdiam=150"
metadata:
  type: project
---

Set up 2026-09-14 for abinitio2D. Source project:
`/mnt/beegfs/elmlund/processed_test_datasets/flipqr/8_extract/flipqr.simple` (owner reboulcf,
group-readable). My work dir: `/mnt/beegfs/elmlund/afan/flipqr`.

**Path trap:** it is `/mnt/beegfs/elmlund/...`, not `emlund`.

**Geometry trap:** the project `log` says `smpd=0.822`, but that is the *movie* pixel size.
Motion correction wrote `_intg` at **smpd=1.3** (xdim 2426 x 2346), so the extracted particles are
**box=180 at smpd 1.3** = 234 A. 2514 micrographs, **689,169 particles**, 84 GB of `.mrcs`.
Read the `stk` segment, never the shell log, for the real smpd.

**Diameter:** measured from `automask2D` on the 100 mini_stream cavgs (per-class radius of the
binarized mask): median diameter 116-120 A, p90 ~165 A (those are doublet/aggregate classes),
equivalent-circular p50 109 A. SIMPLE's own pick-time `ESTIMATED BOX SIZE: 116` px = 151 A agrees.
So **mskdiam=150**. That yields lpstart 15.0 / lpstop 6.8 / cenlp 25.0 via
`mskdiam2lplimits_cluster2D`.

**abinitio2D autoscale:** `SMPD_TARGET=2.67` and `MINBOXSZ=88` in
`src/main/simple_abinitio2D_controller.f90` crop 180 -> **box_crop 88** (234/2.67 = 87.6).
`NSAMPLE_DEFAULT_2D=200000` so the 689k set is sampled, update_frac 0.290, 6 stages.

Only copy `flipqr.simple` into the work dir — the stack paths are absolute, so the 84 GB stays put.
Related: [[batch2-archopt-sigill-login-node]], [[simple-script-yes-compenv-trap]].
