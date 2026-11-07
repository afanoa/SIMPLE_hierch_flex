---
name: empiar-10076-deposited-class-validation
description: "EMPIAR-10076 ground truth: deposited maps at /home/afan/emdb_10076/, 32 components is the measured optimum, classes are effectively 2-D, and the ~2x volume-spread target was wrong"
metadata:
  type: reference
---

Davis et al. deposited classes EMD-8440..8449 sit on the remote box at `/home/afan/emdb_10076/`,
Fourier-cropped to box 128 at 3.275 A. Frame check: our consensus correlates **+0.94** with their
class mean at zero shift, identity orientation, sign -1 (we are protein-dark), so the comparison is
sound and any mismatch is real rather than a convention error.

- **32 components is a measured optimum, not a floor.** Corrected noise: 20 comps 7/10 classes hit,
  mean best match +0.679 discrimination +0.377 at 32; 48 comps is worse on every metric
  (6/10, +0.626, +0.122).
- **The "~2x volume spread" target was wrong** — it was read off the published *figure*, where
  apparent size differences are orientation and rendering. The deposited classes differ by only
  **1.22x** in enclosed volume (1183–1438 kA^3). Volume spread must stop being used as a success
  metric; the deviation-correlation matrix against the deposited maps is the right scorer. (This is
  enclosed volume, NOT the map-triple `spread` in [[flex-pca-downstream-scoring-protocol]].)
- **Their 10 classes are effectively 2-dimensional**: C-group intercorrelates +0.66..+0.98, D-group
  +0.28..+0.93, C vs D -0.56..-0.97, B near-orthogonal. cryoDRGN reports 4 major assembly states
  from an 8-D latent, which agrees. Do not expect 10 independent things.
- **Known defect (2026-08-17):** the 2D shift search on `jclean10076/1_abinitio2D` diverged —
  imported shifts were sane (-37..110 px, 0/799 beyond 128) and came out -444..544 with 35% beyond
  128, so every class average is a blob and the junk-rejection comparison against cryoDRGN is void.
  Rebuild from stack + CTF-only parameters (no angles, no shifts). Trap: `import_particles` writes
  the populated project into `1_import_particles/` while the top-level `.simple` stays an empty
  stub.
