---
name: cls-expansion-fortran-toy
description: end-to-end Fortran cls_expansion toy (subunit motion + periphery + head compositional state) on batch2; 17 s per run; scores amplitude/err per region and head purity/occupancy
metadata:
  type: project
---
scripts/cls_expansion_proto/fortran_toy/ (gen_toy.py, run_toy.sh, score_toy.py; copies in
/mnt/beegfs/elmlund/afan/clsexp_toy/scripts). `sbatch -o OUT.%j.out run_toy.sh OUTDIR [SNR] [N]`; RUNTAG names the run dir,
CLSEXP_ARGS are appended to cls_expansion. Data reused if OUTDIR/toy exists (snr002 = no head, head002 = head 65%).
batch2 python is numpy 1.14 (no default_rng, no fromfile offset=); run_toy needs SIMPLE_QSYS=local after the SIMPLE_* unset.
Baseline 2026-10-07 (build ba2c59520, rank 8, pooled rendering, SNR 0.02): rendered subunit amp 0.90 err 0.13, periphery
0.99/0.05, head 0.92/0.03, head purity 1.000, occupancy cc 0.99. Rank 12/16/24 did not help the subunit (0.86-0.90, err 0.13)
and made the periphery and member averages worse: rank 8 stays.
User requirement (2026-10-07): compositional states (ribosome with/without head) must stay separable - check the head score
after every change. Related: [[python-latent-needs-ctf-weighted-residual]]

2026-10-07 later: kernel weights removed; rendering = placement centres de-shrunk ((I - S)^-1 (zbar_s - zbar)), same scores.
HONEST protocol (subunit_fixes.py / basis_sweep.py: place on noise realisation Y1, estimate from Y2) shows the "seen"
amplitude ~1.0 is selection-inflated; honest subunit capture at lp=10 is only 0.68. Not the linear model (oracle rank-8
basis 0.90), not placement/kernel/residual/local-linear: it is basis estimation noise. Fit band is the lever (6k ptcls,
honest subunit/periphery/head): 10 A 0.68/0.69/0.85, 15 A 0.77/0.83/0.93, 20 A 0.82/0.90/0.96, 25 A 0.86/0.91/0.98,
30 A 0.86/0.91/0.99; 12k at 10 A only 0.78. Fortran toy lp=25 seen 0.91/0.96/0.94 (lp=10 0.90/0.99/0.93); hotspots
(local columns) found 0 at lp>=20. Earlier NOT held-out analysis favoured lp=10 - conflict to resolve on NOT.

2026-10-08: the data-derived band rule (d000ca6d7/9ba0b831a: global = basis even/odd FRC 0.5, hotspots = 0.143) WORKS on
the toy (20 A) but FAILS on NOT: class resolution 2.6-2.8 A (Nyquist), global band 3.7-3.9 A for classes 15/269 (108 s fits)
and 136 A for 181. Suspect: without the pose-tangent nuisance columns, in-plane jitter makes reproducible high-res
"heterogeneity". NOT runs 21 (3 classes, 8.7 min) and 22 (full selection, 37 min) used these bands - treat as invalid.
Spliceosome 10180 classes 49/164/253 (flexcls/10180/clsexp_render3_lp20, lp=20): 1.5 min, hotspots in 164 and 253.

2026-10-09 band growth (user's idea: use the variation FRC to find what reproduces, hotspots from the same curve), scratch
branch band-exp in ~/SIMPLE_clsband (SIMPLE_FLEXCLS_BANDS=grow|frc, SIMPLE_FLEXCLS_SHELLW), local Mac runs. Toys: gen_toy2.py
(sharp blobs sig 0.8 px + 1.5 deg/0.3 px misalignment; TOY_APIX=2.135 for a 4.3 A Nyquist; run_toy_local.sh TOY_MSKDIAM).
Grow = fit at band k, place, two independent half models, pooled FRC of their subclass deviations on every shell, k -> last
shell before FRC < 0.5, repeat (3-4 its from 4 shells); hotspots at 0.143. Deviation error vs truth (same-noise metric):
HR 4.27 A: capped 20 A 0.439, best fixed 15 A 0.345, grow 12.4 A 0.335 (10 A: head purity 0.66 = crowding; grow stops before).
HR 2.13 A: capped 0.420, fixed 6 A 0.315, grow 6.2 A 0.325 (5 A crowds). Soft head002: capped 0.127, grow 0.122.
Per-shell equal weighting HURTS (coarse shells dominate; with grow it collapsed to 136 A). Uncapped FRC rule also fine on
toys - NOT's 5 A first-fit collapse is not reproduced by any toy. Cost: grow 23-56 s/class vs 9-17 s. Not yet checked with
the honest Y1/Y2 protocol (needs the deleted restoration-only mode or an equivalent) nor on NOT class 54.
Growth speed-up (2026-10-09, band-exp): warm starts (each growth step continues from the previous step's extended basis,
joint from joint, each half from its own half; 3 random starts only at step 1), stop when the band would move <= 10%,
never shrink below a band an earlier step confirmed (1-2 random starts collapsed the HR toy to 109-136 A), delivered half
models warm-started from growth's last halves; the final JOINT fit must stay cold (warm final: error 0.335 -> 0.413).
Same bands and scores as full growth on all three toys; class time HR 56 -> 30 s, HR2 55 -> 36 s, soft 23 -> 19 s
(capped 17/9/11 s; most of the rest is the finer final band itself). Env: SIMPLE_FLEXCLS_BANDS=grow GROW_QUICK=yes
GROW_NSTART=3 GROW_WARM=yes GROW_FINALWARM=halves.
NOT 2026-10-09 (32_class_test, classes 193/54/235/95, snapshot bbdde8a0d): reader fix confirmed (read+transform 3-7 s/class,
was 69 s; 4 classes in 2.5-3 min). Growth from 4 shells FAILS on NOT: the 136 A latent only sees intensity/background, its
deviation FRC < 0.5 from shell 2, every class stays at 136 A. From 20 A (GROW_START=27): 54 grows to 6.3 A (FRC >= 0.73 to
Nyquist = reproducible high-res variation, likely residual misalignment, that the FRC cannot tell from motion), 193/95/235
collapse to 136 A because the pooled deviation FRC is LOW at the coarsest shells (0.1-0.4) and higher mid-band. The
"first shell below 0.5" scan is fragile on real curves; the FRC criterion alone does not separate nuisance from motion.
