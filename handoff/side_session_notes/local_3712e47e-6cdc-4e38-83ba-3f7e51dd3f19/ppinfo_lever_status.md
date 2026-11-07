# Per-particle information levers at the arm site (session tree-ppinfo-lever, worktree suspicious-dewdney-bbfa4e, branch claude/suspicious-dewdney-bbfa4e reset onto claude/tree-lever-floor-bccc29 c5176bddc)

Work dir on batch2: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo/ (ppinfo_site.py the numpy port of site_em +
the hold-out acceptance; ppinfo_run.sbatch one config per array task, 4 cores; dumps/<node> = print_project_field text of n.simple and
stageB/b.simple; cache/<node>_B88.npy the particles' Fourier planes at box 88 (1.9 GB per node, 60k particles); runs/<config>_<site>/
summary.json + site_post_<tag>.txt in site_rank's format). Scripts in the worktree under scripts/hetero_tree_proto_ppinfo/, note
doc/implementation_notes/tree_ppinfo_lever_2026_11.md. No Fortran touched.

## 2026-11-02 ~15:20 (batch2 time)
- Read the Fortran (site_em, exec_site_rank, gen_fplane4rec, plane_cell, euler2m, the CTF kernel) and ported it: planes = shifted
  observation x signed CTF on the native lattice of the cropped box, nearest-neighbour class maps with the median-density regulariser,
  cross-half E-steps on FSC-weighted site-masked deviations, top-N hold-out curve over shells kb+2..khi, z at mask_neff.
- Levers as switches: --box/--emlim (band), --maskmode soft|core|sphere|oracle (mask), --poses input|stageB|true (pose),
  --labels em|truth|file (oracle labels = the acceptance ceiling of a perfect classifier; file = rescoring a saved posterior).
- Validation running (job 54871906): the four saved c0_L065s1_site5 posteriors rescored with my hold-out machinery (must give the
  Fortran's site_curve_*.txt: r 0.147 / 0.141 / 0.130 / 0.148 at N* 25000 / 25000 / 15000 / 25000) plus a geometry check (FSC of my
  all-particle NN map against the node's mean/recvol_state01.mrc). Then the EM baseline (must give AUC .75-.78, r .13-.15).
- Facts for whoever follows: pool2 particles are box 160 at 1.464 A (ctf=yes, 300 kV, Cs 2.7, fraca 0.1); site_rank crops to 48
  (smpd_crop 4.88); the true conformer maps are routes/toy55/maps160/docked_all/conf18_emd74209_m160.mrc (rare) and
  conf01_emd74131_m160.mrc (majority); the true poses are detlim/pool2/all.txt (72000 rows, x = y = 0), the fixture input poses are
  refined from a 10-deg / 3-px perturbation of them.

## 2026-11-02 ~15:35 (batch2 time)
- Port validated for geometry: FSC of my all-particle nearest-neighbour map against the node's mean/recvol_state01.mrc 0.96-0.97
  from 117 to 18 A (the shift sign confirmed by flipping it: 0.97 -> 0.56 at 19.5 A). One bug found and fixed by a per-shell
  diagnostic (diag_shells.py): I had stored the signed CTF where the Fortran keeps CTF^2, so the density summed signed CTFs and
  the shells where the ensemble CTF changes sign (15, 20, 23 of 48) came out dead; the plane-slice agreement per ring followed the
  data's CTF lobes exactly (positive to ring 13, zero at 14-15, negative at 16-19), so the planes are right.
- Rerunning the two validation tasks (job 54873350: the saved c0 posteriors rescored + the EM baseline). Phase-2 list (44 configs:
  base/core/stageB on all 7 sites, band steps b64/b64e12/b64e10 on the 6.5 % sites + nulls, truth-label ceilings, true poses,
  oracle mask) is written; it goes out once the baseline reproduces AUC .75-.78 and r .13-.15.

## 2026-11-02 ~16:05 (batch2 time): port validated, ceilings in
- EM baseline reproduces the Fortran (AUC .768/.792 vs .754/.780, r .13-.16, 10 iterations, 133 s per site on 4 cores).
- Ceilings at 6.5 % (L065s1/s2), all at box 48 / 15 A / soft site mask unless said:
  * truth labels (perfect classifier): r .346/.367, z 4.1/4.3 -> the site's map signal PASSES the bar; the loss is the posterior's purity.
  * one M+E step from the truth (perfect references): AUC .880/.883, class>=0.5 1056 at P .81, r .134/.128 (z 1.5/1.4) = the EM's own
    endpoints. At box 48 no classifier crosses with these planes. Needed: top-4000 purity ~0.8 (AUC ~0.98) vs 0.37.
  * same with the EM band to 10 A (box 64): AUC .880 identical -> shells 17-23 carry no per-particle information; hold-out to 7.6 A
    dilutes r to .05 (the truth deviation's FSC crosses 0.5 at ~12.5 A). BAND LEVER CLOSED.
  * core mask: AUC .860 (< .880); its hold-out r peaks at N* 25000 (.33) = nuisance reproducibility of large sets (also z1lo anti-rare
    in the core EM: r .31 at N* 40000, 0 seeds). Null check pending (core_qtoynull/olenull running).
  * TRUE poses + perfect references: AUC .942, class 1862 at P .87, r .405 z 4.9 -> POLE. The pose is a real information lever (the
    fixture's input poses are refined from a 10-deg/3-px perturbation). Realistic version = stage-B poses: EM running (stageB_*),
    ceiling queued (seed1_stageB).
  * whole-molecule E-step (no site mask): AUC .946 -> the site mask under-reads the difference (site 5 = 37 % of the energy); its
    posterior rescored under the site mask queued (rsph). Oracle-mask EM running, oracle ceiling queued.
  * full EM from the truth start drifts to AUC .815, share .12: the EM's fixed point is diluted (the climb no-go seen from above).

## 2026-11-02 ~16:45 (batch2 time): the pose is the lever
- EM with the TRUE poses at 6.5 % (realistic latent-tail starts): z3lo AUC .893/.899, class>=0.5 1743/1960 at P .76/.74, hold-out
  r .297/.316 (z 3.5/3.7 at my n_eff), seed test z 3.27/3.34 -> KEPT POLE at both seeds. Thin margin (1.2x the bar).
- Pose errors of the node's input poses vs truth: median 3.8 deg, mean 9.8, 30 % > 5 deg, 8 % > 10 deg, 4.7 % > 30 deg, shifts
  median 1.6 px at box 160. Stage B: majority 3.4 deg, rare particles 3.75 deg at 6.5 % (no gain) and 3.41 at 8 % (gain: the
  stage-B split is the rare state there) -> stage-B poses are circular as a lever.
- Oracle mask EM: z3lo accepted (r .31, z 3.2 at n_eff 102) but 234/157 seeds at seed z .5/.4 -> no kept pole. Core mask: a
  recalibration (nulls' z rise too), information down. Band: closed (AUC unchanged to 10 A; wider hold-out dilutes).
- Nulls under the seed test: qtoynull base = 0 kept poles (the two small-class endpoints accepted at z 14 have 23-24 seeds at
  seed z ~0, as in the Fortran). 8 %: accepted endpoints, seeds 2588/2940 at seed z 2.86/2.73 -> not kept (Fortran: 2.2-2.3);
  10 %: KEPT (seed z 3.5-3.6).
- Running: OLE null reruns at the right box (70; the first set ran at 48 by my default and had an empty hold-out band),
  the seed-test rescoring of the earlier runs (36), and the pose-accuracy requirement (true poses jittered 1.5 / 3 deg).
- Cost so far ~14 core-h.

## 2026-11-02 ~17:30 (batch2 time): DONE except the commit
- VERDICT: no lever the node can pull crosses the bar at 6.5 % with the nulls pole-free. Band closed by information (AUC unchanged
  to 10 A, wider hold-out dilutes even the truth map). Core mask = recalibration (AUC down, nulls' z up, qtoynull gets a kept pole of
  17 seeds). Oracle mask / whole-molecule E-step raise the information (AUC .927 / .946 with perfect refs) but no kept pole under the
  site mask (z 2.7 at best). POSE is the lever: true poses -> EM kept pole at both 6.5 % seeds (z 3.5-3.7, seed z 3.3, 1.2x the bar);
  requirement ~1 deg median on the rare particles (1-deg jitter z 2.95, just short); node poses 3.8 deg median; stage B helps only
  at 8 % where the state is already in the stage-B states (circular). Budget: perfect refs AUC .88 (d' 1.66) vs needed ~.98 (d' 2.9).
- Nulls: qtoynull and olenull (box 70) clean at baseline, band and stage-B levers; core fails on qtoynull. 8-10 %: purity kept.
- 2-deg random pose jitter makes the EM take reproducible 45/55 partitions (z 3.6-4.1): a hazard for any pose marginal not
  consistent with the images; random jitter is not a model of refined-pose error.
- Note: doc/implementation_notes/tree_ppinfo_lever_2026_11.md; scripts/hetero_tree_proto_ppinfo/ (ppinfo_site.py, ppinfo_run.sbatch,
  make_lists.sh, make_rescore.sh, ppinfo_table.py --md, pose_err.py, diag_shells.py, mrcio.py). Memory ppinfo-lever-verdict-2026-11-02.
