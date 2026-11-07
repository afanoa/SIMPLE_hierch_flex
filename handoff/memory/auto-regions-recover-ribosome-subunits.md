---
name: auto-regions-recover-ribosome-subunits
description: "2026-10-03 EMPIAR-10028 Pf80S: regions_auto (eigengap) picks B=2 and splits SSU/LSU at 91% voxel agreement (SSU recall 0.97, LSU region 98% pure); B=3 splits the SSU in two; our consensus is mirrored vs the deposited models"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-04T01:16:25.622Z
---

**Setup:**
- Script `~/autostates/ribo_subunits.py`, outputs in `10028/10028_simple/RIBO/subunits/`.
- Deposited 3J79 (LSU) and 3J7A (SSU) fitted into `RIBO/ctrl/1_reconstruct3D/recvol_state01.mrc` (box 360 @ 1.34 Å).
- Fit: 12° Euler grid over both hands with FFT translation, then local refinement. Fit cc 0.761 in the automask.
- **Our map is the MIRROR hand** relative to the deposited models.
- Basis: 20-column W_refit (`flex_verify_20260915/W_refit`, polished pcs).

**Results:**
- regions_auto eigengaps: B2 0.081, B3 0.049, B4 0.048, so it picks B=2. Within/between |corr| ratio only 1.17.
- **B=2:**
  - Region 1: 33,867 SSU + 7,358 LSU voxels (82% pure).
  - Region 2: 51,558 LSU + 953 SSU voxels (98% pure).
  - Agreement 0.911; SSU recall 0.973.
- **B=3:** LSU stays one pure region (99%); the SSU splits into two regions (purity 0.56 and 0.73). Agreement 0.786.

**How to apply:**
- Co-variation regions do recover the subunits on a ratcheting ribosome, even at a weak 1.17 ratio.
- The 18% LSU share of the SSU region is the interface/co-moving part; not localised yet.
- Earlier "49/51 split, never validated" ([[toy3-body-frame-test-system]] ribosome note) is now validated.

Related: [[rnasep-goal-is-conformations-not-periphery-resolution]], [[multibody-flex-means-region-pca]].
