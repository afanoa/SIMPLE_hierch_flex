---
name: pose-state-coupling-toy-verdict
description: "2026-10-04 overnight toyPS: coupling poses and states (bin frame correction, +per-state search, +re-embed) does not improve state maps on any toy and worsens latent ordering; consensus pose bias is ~95% coherent per conformation, so state maps are only misframed"
metadata:
  node_type: memory
  type: project
  originSessionId: ba648e52-6783-46a5-8d99-bdeeebdef2ff
  modified: 2026-10-04T06:56:26.128Z
---

User asked (2026-10-04, after "will cluster tracking solve heterogeneity?") to develop the coupled pose/state idea on a
toy overnight and decide if it is useful. Write-up: `doc/implementation_notes/pose_state_coupling_toy_2026_10_04.md`.
Code: batch2 `~/toyPS/` (`loop.sbatch` arm A, `loopB.sbatch` arm B = + per-particle search vs bin maps via
`toyps_em.py --binref`, `oracle.sbatch`, `binscore3.py` frame-invariant scorer, `maps_al.py`, `bias_decomp.py`,
`readout.py`, `imprint.py`, `r3s.sbatch`). Data: `/mnt/beegfs/elmlund/afan/rnasep_regions/toyPS/` (t3null, t3offc, t3lowN,
t32dof; toy3 centred data is toyJ/t3s0.02).

**Findings:**
- **Pose bias is coherent per conformation.** Of refine3D's 4.39 deg median error on toy3, 4.22 is the
  per-conformation mean and 1.42 the scatter (null scatter 1.13).
- **State maps are misframed, not blurred.**
  - Scored after a rigid core alignment, the uncoupled maps are within 0.002-0.006 (8/6 A) of the true-pose
    ceiling on toy3, the off-centre hinge, and 12k particles at SNR 0.01.
  - Frame-free scoring is mandatory for any pose-related map claim. The toyJ "maps 0.75 -> 0.89" was a frame
    artefact.
- **Even the oracle does not help the maps.** Oracle poses (each particle searched against its own true volume)
  improve only the latent ordering: rho 0.942 -> 0.963, and at low N 0.899 -> 0.951; with a curved readout up to 0.985.
- **Practical coupling loses.** Correcting poses by the ESTIMATED conformation turns latent error into pose scatter
  (1.42 -> 1.83 deg), so re-embedding orders worse (0.942 -> 0.91, low N 0.899 -> 0.85). Per-particle search vs
  bin maps drifts further each round.
- **Null is clean for the loops.** No fabrication: bin rotations <= 0.33 deg, pose error flat, selection-bias
  readout flat. But sorting on noise alone gives an extreme-bin difference 2.8x the even/odd noise with
  cross-half 0.45.
- **refine3D_states imprints, and its null states pass the gap test.**
  - On the homogeneous null its 4 labels go Mahal2 0.71 -> 3.15 when re-embedded at its final poses
    (RNase P 0.19 -> 3.74).
  - They pass the cross-half GAP test: +1.28 / +0.97 / +0.55 at 60-20 / 20-12 / 12-8 A, vs +1.32 / +1.17 /
    +0.86 for real states.
  - They are view/defocus splits: NMI(state, 50 view cones) 0.163 vs permutation 0.0008. RNase P's labels:
    view NMI 0.027.
  - So validate refine3D_states states with view/defocus NMI controls, not the gap test alone.
  - On toy3 it assigns states worse than plain latent bins (theta variance explained 0.774 vs 0.841;
    oracle 0.882).
- **Coarse states are the one case where bias blurs maps.** On the 2-DOF toy, bins along one latent direction
  score 0.964 vs truth 0.989. Coupling does not fix it (0.963); 9 k-means 2D states do (0.984).
- **RNase P S4 real ladders show no frame offsets** (<= 0.23 deg).
- **flex_pca null states pass the gap test only in the coarse band:** +1.15 at 60-20 A, but +0.26 / +0.04 at
  20-12 / 12-8 A, against real +1.04 / +0.61.

- **Cluster tracking (2D classes -> pooled 3D latent -> states) is unbiased but blind on toy3** (`~/toyPS/classlink2.py`,
  `cls2d.sbatch`).
  - View-centred, cross-decided class states have the fewest null fakes: fine-band gap +0.02 vs +0.30 per
    particle.
  - Without view-centring, class states are view groups that pass the gap test.
  - But flex2D's 300 heterogeneity-mode classes contain the whole +-15 deg swing at SNR 0.02 (class
    membership explains 7.6 % of theta), so class states recover 3-5 % vs 85 % for per-particle bins.
  - Only compositional states (state-pure flex2D classes) could make it pay; untested.
**How to apply:**
- Do not build pose-state coupling into flex.
- Core-align state maps before reading difference maps.
- Read the flex gap test at 20 A and finer only; add view/defocus NMI controls for refine3D_states states.
- The state-assignment headroom needs conformation estimates better than the latent's own, which is circular.

Related: [[robust-frame-loop-removes-heterogeneity-pose-bias]], [[joint-pose-heterogeneity-toy-verdict]],
[[rnasep-final-poses-make-states-separable]] (refine3D_states circularity test: see the write-up).
