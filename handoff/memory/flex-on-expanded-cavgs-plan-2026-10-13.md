---
name: flex-on-expanded-cavgs-plan-2026-10-13
description: "Plan and results 2026-10-13: flex_pca on cls_expansion averages for discrete RNase P states; flex_view_centre works on S-pure averages, expanded averages are state mixtures, refine=states (hard, fixed pose) freezes partitions, S labels need their own poses"
metadata:
  node_type: memory
  type: project
  originSessionId: 34b66a02-03d2-42df-91f8-915afbf80120
  modified: 2026-10-13T22:33:44.970Z
---

User goal: model the motions at the 3D level with flex_pca on the EXPANDED class averages and get discrete states
directly (no refine3D_states_cavgs march). User decisions 2026-10-13:
- Poses: use solve3D_cavgs single-state on the expanded averages ("I know it works"), not parent poses / docking.
- Noise term: "try both" = (a) least squares with no noise model (drop prior, calibration, XD deconvolution),
  (b) per-average noise scaled by 1/members (sigma_pop-like).
- Implement in Fortran if quicker to test. Go above K=2; "the latent should support quite a few states" -> show the UMAP.
- Expansion may be bumped to 100 sub-classes per class (ncls=100 in the coupled build 9eaef62df counts PER class).
- Work autonomously for 5+ hours; cluster has room; tests must be optimised (distributed arms, parallel).

Runs (batch2, /mnt/beegfs/elmlund/afan/rnasep_coupled/):
- flexcavgs/mkcavgproj.sh NAME EXPDIR ORIS OUTDIR: averages-as-particles project (stktab even+odd stacks, ctf=no,
  per-row pose from solve3D_cavgs final_oris.txt, class = sub-class index, eo 0/1); flexarm.sbatch = one standalone
  flex_pca arm (nparts=10 nthr=8); arms under flexcavgs/arms/<TAG>.
- cavgs16 = full692k_2d 16/class expansion + r3s_cavgs/CTRL_EXP poses/map (1,792 sub-classes, 2,504 active rows).
- full692k_2d100/exp100.sbatch: ncls=100 expansion -> s3d/ solve3D_cavgs (pop>=20) -> cavgs100 -> arms e100_*.
- Step 5 plan: sub-class labels -> particles (cls_expansion_labels.txt) -> selection infile=states (multi-state ok)
  -> reconstruct3D on the 692k particles at the 24_refine3D_auto poses.

Results 2026-10-13 ~18:10 (cavgs16b = cavgs16 with eo alternating by sub-class; flex workers refuse a part holding one eo
half -> "degenerate eo split", GW:99; master then hangs):
- New code (uncommitted ~/SIMPLE_master1013, snap reseed13): params flex_cavgs=no|model|ls. model = per-row contrast
  a=clamp(myth/e_mm) in the embedding (was hard 1 -> latent radius tracked member count, AUC 0.89); ls = also prior 0 +
  1e-3 relative ridge, no MCFA, E[zz']=zz', XD deconvolution skipped (coords = LS z). sigma_pop arms need cls2D in the
  project: cavgs16pop = K8D2 refine3D_states_cavgs work project + replace_project_field ptcl3D/ptcl2D from cavgs16b
  (replace_project only does ptcl2D/ptcl3D/stk with equal counts; no program writes cls2D pop).
- lp20 neigs3 npi2 K8 min_neff=20 (default 2000 prunes averages): current code e/o same-state 0.905, size AUC 0.89;
  model 0.933, radius corr -0.07, size AUC 0.73; ls 0.934. w2_m UMAP = TWO ISLANDS: 75 % / 25 % of particles, e/o same
  island 0.998, view NMI 0.055, 103/112 parents hold both; island-2 sub-classes ~5x fewer members (minority state or size?)
  ls arm gives the same split (96 % row agreement). Particle reconstructions (states2ptcls.sbatch, REMAP old=new:...) pending.

- The two islands are a NOISE-LEVEL artefact: island-2 = small sub-classes, renders 1.8x (12 A) - 3.2x (6 A) more
  high-frequency power, 2x further from the parent mean, same structure (montage); halves share the noise level so e/o
  0.998 means nothing; island maps differ in the core only (cc 0.983 at 12 A), island-2 enriched in junk sinks S3/S5.
- Member-count AUC 0.73-0.83 survives contrast, ls, sigma_pop, XD (now K up to n/100 via min_per_comp, snap reseed14),
  sigma_pop=frc (per-shell 1-FRC noise, frcs=, reseed15), merged averages. Threshold pop>=400 removes it (AUC 0.6) but
  then VIEW takes over (states span 7-61/112 views, NMI 0.33) = the Oct-8 rendering-regime confound.
- Step-5 particle maps (states2ptcls, compare_maps.py vs S1/S2/S4 deviations from their own means): best flex states
  match S1/S2/S4 changes at only r 0.30-0.35.
- Next: parent-referenced input y' = y - ybar_parent + a_p R_i (R from solve3D_cavgs cavgs_reprojs.mrc, interleaved
  avg/reproj pairs), flexcavgs/mkpref.py -> cavgs16pref; wave-8 arms running.

**Why:** refine3D_states_cavgs drains/flips because per-state maps from subsets of the 112 views are
under-determined; flex on averages failed because latent radius tracks member count (AUC 0.86) and view leaks.
**How to apply:** continue from these paths; report UMAPs; see [[rnasep-cavg-states-findings-2026-10-13]].

Results 2026-10-13 ~19:00-19:40 (positive control + view centring):
- S-pure control (flexcavgs/spure, mkspure.sbatch): make_cavgs per 25_refine3D_states S on full692k_2d -> 112 parent
  averages per S (class numbering = original 500-class CTRL_ORIGF numbering, cc >= 0.98), posed with the parent's
  CTRL_ORIGF solve3D_cavgs pose. flex lp20 model: placement = groups of PARENTS (NMI(state,S) 0.00), but the latent
  carries S: parent-centred latent predicts S1 vs S2/S4 at 0.93-0.98 (parent-grouped CV) and k-means on z minus the
  parent mean finds a 98 %-pure S1 cluster (NMI 0.48, e/o 0.90; S2/S4 partly mixed).
- Particle-level flex latent (full692k/1_flex_pca): parent means explain only 2.7 % of the variance; centring does
  not help (S1-vs-rest AUC 0.70 either way). The view offset is an AVERAGES phenomenon (noise averaged away).
- Expanded averages (16 and 100/class): parent-centred k-means NMI(S) <= 0.018 vs the sub-class ceiling 0.035-0.037:
  the expansion sub-classes are S-mixed; centring cannot create state purity the sub-classes lack.
- New code (snap reseed16, uncommitted): flex_view_centre=yes (flex_cavgs only) subtracts the per-group mean latent
  (group = ptcl3D class; averages projects must carry class = parent: mkcavgproj.sh 5th arg CLASSMAP, mkspure CLASSPAR=1)
  before deconvolution/placement, writes flex_pca_umap_vc.*.
- TRAP: cls_expansion_class_map.txt rows are NOT in sub-class order (column 1 = global sub-class id = stack slice =
  cls2D index). readout.py/exp_cluster.py indexed by row until 19:30; parent partition was still right (ARI 1.0, blocks
  of 16/100 align), pop was a rank proxy (log corr 0.84/0.92). Fixed; corrected size AUC 0.85 (16/class), 0.64 (100).
- make_cavgs bug: on projects whose rows are not stack-ordered (full692k_2d) it accumulates only ~40 % of active rows
  (cavger_assemble_sums walks stack row ranges); cls2D pop still reports the full count. Task chip spawned.

Results 2026-10-13 19:15-19:55 (batch2 clock):
- Fortran flex_view_centre=yes (snap reseed16) on S-pure S1/S2/S4: K=3 states, S1 state 216/224 S1 rows, NMI 0.41,
  e/o 0.86 (raw 0.00); 5 S: a 100 %-pure S1 state, S3+S5 together. Centring removes 94-96 % of the latent variance.
- Same on expanded averages (cavgs16b_vc/cavgs100_vc, class = parent): NMI(S) 0.018/0.007, state maps vs S |r| <= 0.36.
  Known S compositions explain 7.6 % of the expanded latent; archetypes land on junk (S3/S5). lp 12/20/30 expansions
  (full692k_2d_lp{12,30}): within-class S1 AUC 0.725/0.728/0.709 -> band is not the limit.
- refine=states (snap reseed17, uncommitted): classification WITHOUT alignment in refine3D (strategy3D_states in
  simple_strategy3D_eval.f90: every state at the stored proj/inpl, argmax; pose/shift never written). Arms
  flexcavgs/rstates/{fv4,rnd4,fv8,S5} (lp 8, nspace 5000, maxits 6) running from 19:55.
- Notes: doc/implementation_notes/planned/flex_pca_class_averages.md in ~/SIMPLE_master1013.
- refine=states first build (reseed17) wrote nothing back: a strategy's srch must call self%oris_assign (the matcher
  calls srch only). reseed18 fixed; reseed19 sets mi_proj = state kept, so multi-state convergence = per-state
  keep fraction > 0.96 (OVERLAP_STATE_JOINT). lp 8 arms (rstates/*_b, converged in 3 its): S5 seed -> 55 % relabelled,
  NMI(final,S) 0.10 (enriched only: 48/61/48 % S1/S2/S4); rnd4 stays random (NMI 0.001); fv4/fv8 mildly enriched.
  => fixed-pose lp-8 per-particle evidence is weak (cf. 2D latent AUC 0.72). lp 5/4 arms (rstates/*_lp5, S5_lp4) next.
- SELF-BIAS: same-half refine=states freezes ANY partition after one reshuffle (random K=4 lp 5: keep 0.27 at it 1,
  0.992 from it 2; S seed: 0.98+ from it 2) - the particle's own half-map carries its noise. reseed20 scores against
  the OTHER half (pftc set_eo flip inside srch_states, restored before sigma2). Cross-half arms rstates/{S5x,S5vx,
  rnd4x,fv4x}_lp5 from ~21:05; S5v_lp5/8 = same-half with the 25-run maps as vol1..5 (iteration-1 keep vs sharp maps).
- Freeze is EINSTEIN FROM NOISE, not self-inclusion: cross-half random K=4 still keeps 0.993 from it 2 (both halves
  sorted against the same it-1 refs). Keep ~1 never validates a hard partition; use seed/half reproducibility.
- S labels need their OWN poses: vs the 25 maps at consensus poses keep S1/S2/S4 ~0.55-0.64 at lp 5/8/15/20 (band-
  independent); at the 25 run's per-state poses 0.93-0.97. Consensus vs 25 poses: median 5.4 deg (S3/S5 35-54 deg).
  => NMI vs S under-rates every fixed-pose method; part of S is the joint pose-state fit. Next: soft refine=states
  with a local per-state pose window. Note: doc/implementation_notes/planned/flex_pca_class_averages.md (master1013).
