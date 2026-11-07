---
name: null-control-fixtures-2026-10-27
description: "No-heterogeneity controls (user 2026-10-27): qtoy2-matched null (routes/toy55/simqtoynull, one conformer) and OLE-matched null (null_controls/ole_null_v1, consensus reprojected at OLE poses/CTF, measured noise); every tree/rule change must not split them"
metadata:
  type: project
---

User 2026-10-27: "we should implement a control where there supposed to be no heterogeneity."
- qtoy2-matched null: `routes/toy55/run_qtoy_sim_jf.sbatch` (copy of run_qtoy_sim with JF env) QT=qtoynull, CSV qtoynull.csv = conf01
  only, JF=0.0002 (12 blob + 12 broken token junk), N 60000, SNR0 0.0012, same B/q quality model and real RNase P poses/CTF as
  qtoy2; refine3D from 10 deg / 3 px perturbed truth -> simqtoynull/src.simple, truth simqtoynull/truth_qtoy.txt (job 54385675).
- OLE-matched null: `hetero_tree_regression_20261026/null_controls/ole_null_v1/` (ole_null.py prep|measure|merge + ole_null.sbatch,
  job 54385840): 5_refine3D_auto/rec_final_state01.mrc Fourier-cropped to box 180 / 1.908 A, simulate_particles noise-free at
  OLE's refined poses + CTF (80,599 rows), noise and amplitude spread (qsig) MEASURED on 2,000 real OLE particles cropped to the
  same frame (measure.json), B ~ U(0,200); refine3D c2 from 10 deg / 2 px perturbed truth -> src.simple.
- Tree launcher datasets `qtoynull` (minpop 900, maxdepth 7) and `olenull` (c2, mskdiam 240, minpop 3000, maxdepth 20) in
  amp_mstep_v1/tree.sbatch. Ideal: one genuine state; any accepted split is a false positive (check snrq/bfac association).
Related: [[efficient-testing-over-wide-queues-2026-10-27]], [[flex-pca-amplitude-mstep-omission-2026-10-27]].
