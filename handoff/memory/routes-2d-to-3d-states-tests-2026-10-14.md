---
name: routes-2d-to-3d-states-tests-2026-10-14
description: "2026-10-14 evening tests of ways to turn cls_expansion 2D sub-class heterogeneity into 3D states: pose-binning refuted on RNase P (1 reproducible direction per view), label-proportion (Sinkhorn) E-step recovers a missing state with oracle bag proportions but needs |err| <~ 0.06, per-view two-endmember unmixing + 3D-consistency sync gives a map-defined sub-class coordinate with |corr| 0.74 to the S1 share on RNase P; pv3 kernel toy splits on pose (no pose mean)"
metadata:
  type: project
---

Scripts: `scripts/cls_expansion_proto/routes/{rnasep_tangent_hysime.py, view_count.py, unmix_sync_pv.py, unmix_sync_dump.py,
pv_bag_estep.py, run_*.sbatch}` and `scripts/cls_expansion_proto/subclass_pool/bag_estep.py`; batch2 work dirs
`/mnt/beegfs/elmlund/afan/routes/{tangent,viewcount,unmix_pv3,unmix_rnasep,pvbag}`; local toySE logs `~/toyse_snr0{3,1}/bag_estep*.log`.

- **Pose binning refuted (RNase P avg_dump, 77 classes, 16/class pose-mean expansion):** 6 pose-tangent columns explain
  0.10-0.15 of the sub-class deviations vs 0.09-0.10 for tangents at a random view; the stored member poses (pm) 0.06;
  even/odd reproducibility 0.77 / 0.86 / 0.26 / 0.07 at 117-30 / 30-20 / 20-12 / 12-8 A barely changes after pose removal.
  ONE reproducible direction per view in every band (69-71 of 77 classes); kernel100 (100/class) shows a second direction
  at 30-20 A in ~40 % of classes but kernel siblings share members (even/odd correlated), so that count is inflated.
- **Label-proportion (Sinkhorn per bag) E-step, toySE SNR 0.03 / 0.01, bags = pose-mean sub-classes, seed lacking state 13:**
  plain E-step 0.79 / 0.68 (13 never found); ORACLE bag proportions 0.99 / 0.89 (13 recall 0.99 / 0.81, recovered from
  nothing); proportions with |err| 0.06 (Dirichlet noise) 0.86 / 0.79; shrunk 30 % to global 0.82; shrunk 60 % 0.68;
  bag means of the first E-step (|err| 0.15) 0.71; archetypes on the pooled 3D latent (|err| 0.22) 0.49. The lever is real
  but needs per-sub-class proportions accurate to ~0.06; the 3D latent cannot supply them. On pv3 particles (sub-classes
  92 % pure, strong per-particle evidence) the plain E-step from the sub-class-majority seed already reaches 0.970 and the
  oracle-proportion Sinkhorn E-step 0.968: the constraint only pays when the seed is poor and per-particle evidence weak.
- **Per-view unmixing + synchronisation (route 2):** K must equal (reproducible directions + 1) = 2 per view; K=3 on
  three states along one axis gives negative R2. RNase P (dump, lp 15): the per-view first direction correlates with the
  S1 share at median |corr| 0.65 (50/77 views > 0.5); near-view edge sign agreement 0.72; spectral Z2 sync gives a global
  coordinate at corr 0.40 (local-sign agreement 0.85); 3D consistency (two endmember maps reconstructed from all views,
  reproject, Hungarian per view) synchronises from a RANDOM start to 0.30 and links the 4 view-graph components. The
  MAP-DEFINED axis (each sub-class average projected on P_B - P_A at its own pose) reaches |corr| 0.74 with the S1 share
  (Spearman 0.68), S1 x1.6 in the S1 fifth; ANY re-reconstruction from the coordinate degrades it (soft weights
  0.74 -> 0.56, view-balanced hard extremes 0.74 -> 0.43): keep the first endmember maps. Bands: lp 15 0.74, lp 20 0.68,
  lp 30 0.70; all 86 views (min 50 members, spread 20 deg) 0.74. Endmember map difference vs S1 - (S2+S4)/2 per band
  117-30 / 30-20 / 20-15 A: 0.58 / 0.36 / 0.12 (lp 15), 0.62 / 0.48 / 0.07 (lp 20, all views); particle-latent unmix
  route was ~0.70 at 15 A. Calibrated to the S1 share these proportions would carry |err| ~0.15: below the 0.06 the
  Sinkhorn E-step needs, so the 2D coordinate is a seed/ordering, not yet a proportion supervision.
- **pv3 kernel toy (no pose mean, 4-deg view clusters):** the leading reproducible direction per view is the sub-class
  POSE offset (|corr| 0.55-0.69 with de1/de2/de3/dx/dy, 0.21 with the truth axis); removing the consensus at each
  sub-class's own pose leaves |corr| 0.18 with truth; the route cannot be judged on this toy without a pose-mean expansion.

- **pv3 redone with the POSE-MEAN expansion (coupled build 9eaef62df, `cls_expansion vol1=... ncls=16` + SIMPLE_FLEXCLS_DIRECT_AVGS,
  chain `scripts/cls_expansion_proto/routes/run_pv3pm.sbatch`, dirs routes/pv3pm/pv3pm_{true,map}):** the pose mean fixes the
  2D partition (sub-classes 99 % pure, ceiling 0.993 with true poses / 0.985 with map2ptcls poses, vs 0.902 kernel without
  it; enrichment x90). Route 2 still fails there: the truth-predicted deviation (truth volumes moved into the solve3D frame
  through the docking transform) carries only 3-9 % of the reproducible within-view sub-class variation, the rest is
  pose/alignment/CTF nuisance that reproduces between halves; nuisance regression (mean per band + tangents), a 50-member
  floor, per-shell whitening, an envelope-invariant ratio coordinate and size-EQUALISED re-averaging (25+25 members, one
  regulariser, `equal_avgs_pv.py`) all leave the per-view first direction at |corr| 0.2-0.3 with the state. An oracle
  axis (maps from truth-extreme sub-classes) reads the state per view at 0.5-0.6 but the axis EM is not self-consistent
  (drifts to 0.1) and random starts never find it. Subtle 20-deg motion: the state is a minor part of the per-view
  scatter; RNase P (S1 vs S2/S4, per-view |corr| 0.65) is the regime where the per-view step works.

- **Joint 3D fit (flex on averages) on the pose-mean pv3 expansion (`run_pv3pm_flex.sbatch`, `avg_latent_cluster.py`):**
  on the delivered DIRECT averages flex's gmm placement collapses to one state (particle accuracy 0.50 all sub-classes,
  0.55 with a 50-member floor, 0.55 on size-equalised 25+25 averages, 0.57 at K=2) and the cross-half E-steps drift down
  (0.55 -> 0.52). The latent itself carries the states only moderately there: with the floor, LDA CV 0.90-0.92, k-means
  0.78-0.81 (sub-class level), composition ridge R2 0.44-0.67; on the kernel pv3 reference the same readout gives LDA
  0.991, k-means 0.975-0.984, gmm placement 0.988 (particles 0.897 = ceiling). So the kernel averages (overlapping kernels,
  neff ~280) are the better input for the joint fit; the coupled build's direct hard-partition averages (50-260 members)
  are noisier, and flex's gmm placement fails on them where k-means still gets 0.81. A pose-mean expansion WITH kernel
  averaging does not exist in one build (pose mean = coupled 9eaef62df, kernel = master 5e1e6aeee). k-means labels on
  that latent (`kmeans_sublab.py`, SUBLAB= in the chain) carried to the particles give 0.82 (state 1 recall 0.32), map
  differences vs truth 0.55-0.59, and the cross-half E-steps DRIFT DOWN 0.82 -> 0.74 (fractions move to the minority,
  purity falls); the kernel pv3 seed (0.897) rose to 0.95 under the same E-steps. Toy verdict: on pv3 the kernel
  expansion + flex on averages + E-step is the working chain; the pose mean buys purity but its direct averages lose it.

- **RNase P 100/class (2026-10-15, `kernel_to_dump.py` adapter + `run_*_route2.sbatch`):** the KERNEL100 expansion (master
  build, no pose mean) is dead for the route: per-view S1 |corr| median 0.15, maps nothing, in every variant (spread 30,
  pose-residual, lp 15/20); tangents explain 0.17 vs 0.09 null, not full pose binning, the leading direction is just not
  the state. The POSE-MEAN 100/class expansion (`rnasep_coupled/full692k_2d100`, coupled build, direct averages, pops
  p10/50/90 6/24/124) works like 16/class per view (S1 |corr| median 0.60-0.65 over 69-83 views) but the synchronised
  map-axis coordinate is a bit lower (0.68 vs 0.74) and the map difference vs S1-(S2+S4)/2 0.52/0.29 vs 0.58/0.36; a
  third of its views show 2-3 reproducible directions at 117-20 A (the K=2 model ignores them). Figures
  `routes/unmix_rnasep/fig_endmembers_{16,100}perclass.png` (projections of the endmember maps vs the S1 reference,
  in-band corr 0.49 / 0.43): the large S1 rearrangement appears as a block-shift dipole in 3D at 117-20 A. Verdict: the
  pose mean is what makes the per-view step work on RNase P; more sub-classes do not sharpen the large motion.

- **Coordinate ends reconstructed from PARTICLES (`run_particle_ends.sbatch`, 2026-10-15):** view-balanced groups (every
  class gives the same fifth of its particles, by whole sub-classes, from the low / middle / high end of the 16/class
  map-axis coordinate; `axis_coord_first.npz` = the coordinate off the FIRST endmember maps, the only one to carry
  forward): low 167k particles S1 share 0.39, high 169k S1 0.19, mid 236k 0.28; each map 5.0-5.4 A at FSC 0.143 (box
  160, 24_refine3D_auto poses); the high-low difference reproduces between half-sets 0.57 / 0.39 / 0.03 at 117-30 /
  30-20 / 20-12 A (difference-map FSC 0.143 = 19.5 A) and matches S1-(S2+S4)/2 at 0.70 / 0.57 / 0.18. So the 2D-seeded
  groups give 5 A maps whose DIFFERENCE is defined only to ~20 A: the ends are ~40/60 mixtures (coordinate-S1 corr 0.74
  at the sub-class level, sub-classes themselves mixed), and mixtures average the fine change away. Maps downloaded to
  ~/Downloads/rnasep_route2_endmembers_20261015/particle_ends_16perclass_20pct. Purer extremes (5 % per class: low 68k
  particles S1 share 0.48, high 87k S1 0.19, maps 6.5 A) do NOT sharpen it: difference reproducibility 0.53 / 0.41 /
  0.09, difference FSC 18.0 A, vs S1 ref 0.64 / 0.55 / 0.15. Sub-class-level grouping saturates at ~20 A; the remaining
  lever is a 2-state refine3D seeded from these ends (per-particle E-steps + per-state poses), not run yet (hours).

**How to apply:** for 2D-derived states use K = ndir + 1 endmembers per view, synchronise by 3D consistency, read the
coordinate off the map axis (hard view-balanced extremes, never soft weights); feed per-sub-class proportions to a
Sinkhorn-constrained E-step only if they are accurate to ~0.06. Related: [[averages-unmixing-overnight-2026-10-14]],
[[rnasep-sync-od-coordinate-2026-10-14]], [[view-imbalance-inflates-between-state-metric]].
