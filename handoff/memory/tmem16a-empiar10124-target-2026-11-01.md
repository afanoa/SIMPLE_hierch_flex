---
name: tmem16a-empiar10124-target-2026-11-01
description: "user chose TMEM16A LMNG (EMPIAR-10124) as the \"current methods fail, ours works\" membrane target; download in beegfs EMPIAR_downloads/10124"
metadata:
  node_type: memory
  type: project
  originSessionId: 851a7b77-9a45-48b1-b2b6-396c1a7b5829
  modified: 2026-11-01T05:22:11.145Z
---

User (2026-11-01) picked mouse TMEM16A in LMNG, EMPIAR-10124 (Dang et al. Nature 2017), as the heterogeneous membrane target where continuous latent methods have not been applied. Asked me to monitor the download and process it autonomously.

- Data: /mnt/beegfs/elmlund/EMPIAR_downloads/10124/data/16A_LMNG/*.mrcs (box 256, 1.02 A/px, normalized, dose-weighted), plus LMNG_all.star / LMNG_final.star; a Fab-bound set (LMNG_Fab_*) left out at first.
- Conditions (EMD-7096): 1 mM CaCl2 (Ca-bound), C2 imposed, 3.8 A from 251,851 particles (LMNG + Fab-subtracted combined), RELION.
- Hypothesis: per-subunit TM6 bending -> asymmetric dimers blurred by C2 consensus; compare with cryoDRGN on published C2 poses as the baseline.
- Sister set EMPIAR-10123 (nanodisc, EMD-7095) for an environment comparison later. MsbA 10799 and gamma-secretase 10194 have no particles (rejected for effort).

**Why:** user wants a molecule where current methods fail and the tree + flex_pca works.
**How to apply:** process without imposing symmetry (pgrp=c1), on the cluster; check STAR columns before import.
