---
name: rnasep-thirds-selection-bias-null
description: "2026-10-03 S4: refine3D_auto thirds along a RANDOM spectrum-matched direction differ 0.132 in region 1 vs 0.20 real, 0.0003 random thirds; coarse <15 A ladder SNR is at bias level; the real part is the coherent arm change"
metadata:
  node_type: memory
  type: project
  originSessionId: 225cbe63-9384-4d98-963e-b85c64a69e1a
  modified: 2026-10-03T23:54:26.575Z
---

**Control:**
- Thirds along null direction 1 of region 1 (`rnasep_autostates/S4/nullthirds/r1_t{1,3}`), refined with refine3D_auto from the consensus start (`~/autostates/as_nullthirds.sbatch`).
- Scored by `~/autostates/as_r1score_null.sbatch` with the G10 region masks.

**Result:** region-1 1−corr at 15/8 Å.

| comparison | 15 Å | 8 Å |
|---|---|---|
| null thirds | 0.132 | 0.132 |
| real thirds | 0.20 | 0.20 |
| random thirds | 0.0003 | 0.0026 |

The null thirds refine to 3.55 Å, the same as the real thirds.

**Interpretation:**
- About 2/3 of the 1−corr is selection bias: particles sorted by their own images carry their noise into the reconstruction.
- On renders the real difference is one coherent arm change; the bias is patchy speckle. 1−corr can't tell them apart, so use null-matched metrics or shape.

**Ladder level (compose route, 3 null directions per region):**
- Below 15 Å, real axes reach SNR 320–530 against null 310–395.
- In 15–8 Å, real axes reach 7–10 against null 4.3–5.2, i.e. z 5–9.
- Posterior coordinate variance: null 68–102 vs real 69–139, so the latent is mostly noise.

**How to apply:**
- Any state claim needs this null.
- The proposed fix is to sort on one band and evaluate in a disjoint band: noise is independent across shells, so there is no bias.
- That needs a lower band limit in `covariance_kfromto` (src/main/flex/simple_flex_pca_em_fit.f90:496, hard-coded 1) and the E-step ring range. The user wanted that done attended.

Related: [[flex-state-maps-selection-bias]], [[rnasep-multibody-region1-split-is-real]], [[rnasep-goal-is-conformations-not-periphery-resolution]].
