---
name: synthetic-toy-traps-noise-poses-2026-10-29
description: "TRAPS building low-SNR toys: white noise at the measured total variance has the wrong spectrum (calibrate per band); refine3D cannot hold poses on PfCRT-like SNR (10 deg / 2 px start -> 41-115 deg); refine=neigh scans 100 % of the space; refine=cont needs eo + proj > 0; project names cannot contain '.'"
metadata:
  node_type: memory
  type: project
  originSessionId: c82ee840-ba2c-4b73-a246-e79b0fb2c726
  modified: 2026-10-30T00:34:44.652Z
---

2026-10-29, PfCRT micelle fixture ([[pfcrt-micelle-fixture-tree-pose-scatter-2026-10-29]]):
- **Noise spectrum.** White noise at the measured total variance (OLE-null recipe) gave band correlations 0.22 / 0.14 / 0.05 / 0.03 vs the real 0.15 / 0.10 / 0.07 / 0.07.
  - Fix: per-band noise power from the real band correlations (`micelle_fixture.py calib` / `merge_cn`). It matches to about 0.005.
- **Poses.** From 10 deg / 2 px perturbed truth, refine3D made poses worse: prob 41-47 deg (white noise), 113-115 deg (calibrated); neigh with a protein mask 106 deg.
  - This is not a convention bug: true poses reconstruct to Nyquist.
  - The 10 deg / 2 px start alone limits the map to 16 A (the shift error dominates).
  - Calibrate the input pose error to the real input's FSC instead (here 3 deg / 0.75 px).
- **refine=neigh** in recipe17d reported "% SEARCH SPACE SCANNED 100": it is not a local search.
- **refine=cont / fixposes** refuses imported poses: it needs eo in {0, 1} and proj > 0. Import does not keep the state column ([[import-particles-resets-state-and-eo]]); set labels with `prg=selection infile=`.
- **new_project** rejects punctuation in projname ('.' in `cal_a2_s0.5`) and needs SIMPLE_QSYS in the job env.

**Why:** each cost a relaunch; the first two silently invalidate a toy.
**How to apply:** check the realism band correlations and the input consensus FSC before any tree run on a new toy.
