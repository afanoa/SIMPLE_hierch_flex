# IgM: our workflow on DJ's .workspace_2 extraction (pentamer vs hexamer)

Working dir: `/mnt/beegfs/elmlund/afan/igm_cryolo/dj_ws2/` on batch.ncifcrf.gov. Build `~/SIMPLE_master1030`. DJ dirs read-only.

## 2026-11-05 12:00 — input and 2D launch
- Input: DJ `.workspace_2/6_extract/workspace.simple` (NICE keeps one project per step; the top-level workspace.simple is 1 kB).
  - 2,251,109 particles in 34,422 stacks, box 600, 1.09 A/px, absolute stack paths, ptcl2D/3D untouched (class 0, state 1). CTF from `.workspace_2/3_ctf_estimate`.
  - Copied to `dj_ws2/dj6_extract.simple` (pristine); solve2D runs on the copy `dj_s2d.simple`.
- solve2D ncls=300 mskdiam=560, 20 workers x 32 threads, norm: job 54955900 (`s2d.sbatch`, run dir `1_solve2D`).
- 2D fold script `cavg_fold.py`: per class, re-centre on the ring (minimise n=1 angular power share at 90-170 A), report n5/n6 share; full-res labelled montages, 100 classes per PNG. Test on our cavgs_iter031: job 54955931.
- Fold script validated on our cavgs_iter031: tag rule = share n5 (or n6) > 0.25, > 1.5x the other, and ring-harmonic amplitude / image SD > 0.2. Tags only 92 (n5 .50/n6 .13), 138 (.40/.17), 246 (.41/.25), 20 (.28/.08): all 5-fold. No 6-fold-dominant class in our 2D (so our 6-fold 3D is not from hexamer 2D classes). Modulation alone flags ice/carbon gradients (centre offsets 20-180 A), hence the combined rule.
- `post2d.sbatch [cavgs]`: cls2D pops, ptcl_class_state.txt, fold table + montages in `fold/`. Watcher on solve2D: fires at iteration 16 or job end.
- 13:10 solve2D at iteration 13 after 1 h (avg class res 13.2 A, 200k sampled per iteration). Fold on iteration 12 cavgs (2.725 A/px, job 54957362): NO star class yet, no 5f/6f tag; classes = round fuzzy blobs, side-view rods, many ice/carbon gradients. Highest n6 share 0.28 (class 293) at modulation 0.04 = noise.

## 15:50 — 2D result: no star classes, so no pentamer/hexamer split in 2D
- solve2D COMPLETED (job 54955900, 3 h 38 min, 2,251,109 particles all active, final cavgs at 600 px). post2d 54960548: fold table `fold/dj_cavgs_iter031_fold.txt`, montages `fold/dj_cavgs_iter031_p1-3.png`.
- No class passes the 5f/6f tag. Best ring modulation 0.02-0.05 (ours class 92: 0.46). Highest shares: 5-leaning 171/199/111 (0.22-0.28), 6-leaning 124/272/132 (0.20-0.23), all at noise modulation.
- Zoom (cavg_zoom.py, job 54960815; `fold/dj_zoom.png`): ~1/3 of DJ classes are featureless ~250 A discs (raw and 15 A low-pass), next to our 92/246 that are crisp 5-hole stars. Rest: side-view rods/tilted views, ice/carbon/gradient junk.
- Mac: ~/Downloads/igm_DJws2_solve2D_300cls_iter031_p1-3.png, igm_DJws2_disc_classes_vs_our_stars_zoom.png. Waiting on the user for how to proceed (class choice / sub-classify discs / 3D).

## 16:00 — disc re-classification + pick-overlap check (user: run both)
- Disc set = 74 shape-rule classes (edge sharpness < 0.055, aspect < 1.25, flat background; discs 0.041-0.050, rods/ice/gradients >= 0.06; `fold/dj_shape.txt`) + 71, 86, 227 = 77 classes, 611,985 particles (`classes_disc77.txt`, `flags_disc77.txt`). Left out 81 (lumpy), 95, 282 (ice balls).
- disc2d.sbatch: copy dj_s2d.simple -> dj_disc.simple, cls2D selection, count gate 611,985, solve2D ncls=100 mskdiam=560 nparts=10 (job 54961009, run dir `2_disc2D`).
- overlap.sbatch (job 54961013, `overlap/`): pick_overlap.py, DJ vs our template picks per micrograph, 50 A match radius, per DJ group (disc77/other) and per our group (star3 92/138/246, core15, keep175, drop110), both directions.

## 16:20 — pick-overlap result (jobs 54961013, 54961060; `overlap/`)
- DJ picked 34,422 of our 43,331 micrographs (median 57 picks/mic vs our 47). No systematic coordinate offset between pick sets (median -1,-1 px); matched picks |offset| median 21 A.
- 53 % of DJ picks have one of ours within 50 A; 22 % are > 200 A from any of ours (different objects).
- On shared micrographs, share of OUR particles with a DJ pick within 50 A: star3 (92/138/246) 34 % (7,363 of 21,378); core15 62 %; our disc classes 51 %; other kept 62 %; dropped junk 71 %. => DJ's picker under-samples the face-on star views ~2x relative to everything else and over-samples junk. The 7,363 star particles DJ did pick scatter over many DJ classes (largest 152), 39 % in DJ discs -> no star class can form.
- Discs are not DJ-specific: same shape rule on OUR 2D = 51 disc classes, 397,342 particles (16.6 % of 2.40M; DJ 26-27 %). 138,397 (20 %) of OUR final 683k consensus set came from our disc classes (all 'keep'). Possibly relevant to the pseudo-hexamer consensus (not tested here).
- DJ disc particles: 62 % have no pick of ours within 50 A (DJ other classes: 41 %). Matched ones: 46 % in our disc classes, 43 % in our other kept classes, 2.3 % core, 1.2 % star.
- Correction to cls_all_pops.txt: pops there sum to 2x the particles (4.80M); use per-particle class files for counts.
- Disc re-classification 54961009 running (iter 3); watcher b86kf4odz submits postdisc.sbatch (fold + montages in fold/disc_*) at its end.

## 16:45 — partition slip (caught by the coordinating session)
- The disc 2D workers (10 x 32 cpus, `2_disc2D`) run on csbdevel, not norm: disc2d.sbatch's update_project passed only job_memory_per_task, and update_compenv then takes `$SIMPLE_QSYS_PARTITION` (= csbdevel, exported by ~/.bashrc on batch) instead of keeping the project's norm. Left running on the coordinator's say-so (csbdevel idle).
- dj_disc.simple now carries csbdevel in compenv. From here on every update_project passes `qsys_name=slurm qsys_partition=norm job_memory_per_task nthr nparts`. postdisc.sbatch has no update_project and runs on norm itself.

## 18:40 — disc re-classification: no star class
- solve2D ncls=100 on the 77 disc classes (611,985 particles): job 54961009 COMPLETED in 2 h 37 min (31 iterations; workers on csbdevel, see 16:45). postdisc 54964070: `fold/disc_cavgs_iter031_fold.txt`, montage `fold/disc_cavgs_iter031_p1.png`.
- No class passes the 5f/6f tag.
  - Highest modulation: 56, 9, 59, 69, 55, 28 (mod5 0.28-0.34, n=3 dominant). Each is a dense knotted object beside a faint disc; the off-centre bright object makes the harmonics, not a star.
  - Largest n5 share: 0.23 (class 8, blob plus dense dot, mod 0.14). Largest n6 share: 0.23 (class 91, featureless blob, mod 0.02).
- By eye (central 400 A, 8 A low-pass, stretched on r < 170 A):
  - featureless mottled discs, 200-250 A (most classes);
  - smooth featureless ellipses (3, 7, 13, 60, 64, 68, 78, 79, 86, 91);
  - small dense rods, 200-250 x 40-50 A (e.g. 1, 2, 4, 14, 20, 23, 24, 33, 35, 52, 53, 58, 74, 75, 84, 88, 97), possibly Fc-core edge-on side views;
  - faint disc beside a dense object (6, 8, 9, 12, 15, 25, 28, 34, 42, 44, 46, 55, 56, 59, 67, 69);
  - ice/carbon gradients (30, 31, 36, 54, 61, 62, 80, 94).
- Mac ~/Downloads:
  - igm_DJws2_disc_reclass_100cls_cavgs_iter031.mrc (class order = the table);
  - igm_DJws2_disc_reclass_100cls_inner400A_8Alp.png (judging view);
  - igm_DJws2_disc_reclass_100cls_montage.png (600 A; cores saturate).
- Waiting on the user's core/kept/dropped choice. 3D: solve3D with nstates > 1 and the full mask (automsk=no), since a single state drifts to the 6-fold register average (consensus finding, 2026-11-05). Every update_project passes qsys_partition=norm.

## 19:30 — user: "just proceed autonomously"; class pick and 3D launch
- Reference: our core-15 rendered the same way (440 A crop, 8 A low-pass, inner-disc stretch). IgM side views are short 2-5-bead rods (Fc core edge-on, 130-200 A; Fabs averaged out); stars 92/138/246 fill the tile.
- DJ 300-class pick (my call, by eye against that reference; `classes_dj_core51.txt`, `classes_dj_kept48.txt`, flags `flags_dj_core.txt`, `flags_dj_super.txt`):
  - core: 51 clear bead-rod side views, 409,710 particles: 6 7 14 15 22 25 26 28 30 38 39 48 55 57 59 66 84 87 89 93 104 116 117 119 123 125 137 143 157 161 162 165 166 174 194 207 209 213 230 232 240 243 259 260 264 267 275 286 287 289 298
  - kept: 48 other IgM-like classes (tilted/knotted views, rods beside diffuse density, mottled discs with a central knot), 291,777 particles: 3 5 29 31 36 51 52 56 61 64 65 76 81 83 85 107 108 111 114 121 139 140 149 152 156 159 163 164 178 190 192 196 215 218 234 236 244 247 250 252 262 263 270 279 285 288 294 296
  - dropped:
    - the 77 disc classes, 611,985 particles;
    - 124 classes of smooth blobs, ice, carbon edges, gradients, and dot-plus-diffuse objects, 937,637 particles.
  - So ~31 % of DJ's 2.25M picks look like IgM.
- 3D job 54965177 (`dj3d.sbatch`, run dir `3D_dj/`):
  - copies of dj_s2d.simple: dj_core.simple (51 classes) and dj_super.simple (99 classes); count gates 409,710 / 701,487;
  - update_project with qsys_partition=norm, job_memory_per_task 48000, nthr 24, nparts 10;
  - solve3D pgrp=c1 mskdiam=560 automsk=no nstates=3 on the core; master 96 GB;
  - symfold_rc.py per state map (tree/symfold.py plus ring re-centring);
  - solve3D_addon projfile=dj_super projfile_frozen=<solved core> addon_diag=yes, then a fold check on the addon maps.
- Fold-script check on known maps (tb3d_s1/s3, EMD-13921): job 54965178.

## 19:50 — the 5/6-fold reading depends on centring (method question sent to the IgM session)
- symfold_rc.py check (job 54965178) on maps with known readings fails. It reads EMD-13921 ring-centred as n5 .23 / n6 .54 at 100 A, and tb3d_s1 as n5 .47.
- Re-reading the saved projections (tree/igm_symfold*.npz) with five centrings: image, local centre of mass within 60 or 90 A, smoothed peak, min-n1. None reproduces the session's EMD-13921 n5 .60/.76 at 90/110 A.
  - EMD-13921 at 110 A reads n6 >= n5 in all centrings except the image centre (.12/.15).
  - The IgM Fc pentamer is a hexagon with a J-chain gap, so n5 vs n6 shifts with the centre a few A either way.
- Sent one question to the IgM session (local_126e2179) about its "local centring". Until it answers, DJ's state maps get tree/symfold.py (COM, as the brief says) plus projection images to judge by eye.
- 20:05 solve3D workers (10 x 24 cpus) confirmed on norm; stage 1 running. Watcher b5v7m1hci fires at SOLVE3D_DONE / job end.

## 20:15 — coordinator go-ahead: random-300k replicate on csbdevel (job 54965797)
- The IgM session (local_126e2179) asked for the reps/rep.sbatch protocol on DJ's particles:
  - seeded random 300k of all 2,251,109 particles, no 2D preselection;
  - box 200 (3.27 A/px);
  - solve3D c1 mskdiam=560 automsk=no nstates=6 lpstop=10 nparts=10 nthr=24, workers on csbdevel;
  - report population, ring-centred n5/n6 at 110 A and the J-chain Fab rod per class, to compare with rep_A / rep_B.
- `dj_r300k.sbatch` (master on csbdevel, 16 cpus, 64 GB):
  1. Copy dj6_extract.simple to `pool200/dj_pool.simple`.
  2. scale_project newbox=200 for all 2.25M into `pool200/stacks_sc`, ~360 GB. It runs as 16 local parts on the master node, so no worker time budget and no flood of the csbdevel queue.
  3. Seeded draw `seq 2251109 | shuf --random-source=<(yes) | head -300000`; selection + prune_project on a copy in `pool200/r300k/`; count gate 300,000.
  4. update_project slurm/csbdevel, then solve3D (6 states); populations; tree/symfold.py into `r300k/symfold_dj.npz`.
- Order differs from the coordinator's text (prune then scale): scaling the full set first gives the same 300k at box 200, mirrors how rep_A/B came from pool200, avoids a 432 GB box-600 pruned copy that could never be deleted, and leaves a reusable DJ pool200.
- The class-based 3-state job 54965177 (norm, 256 cpus) keeps running unless the coordinator wants it stopped.
- 20:30 the IgM session's method arrived, with a correction.
  - Its ring readout centres by MAXIMISING n5+n6 at 70/90/110 A (audit/fold_rc.py). The ring number depends on centring, so it can't settle 5 vs 6 alone.
  - Centring-free test (audit/fold_ft.py): angular power of the 2D Fourier amplitude |F| of the disc-normal projection. A 5-fold object gives n10, a 6-fold one n6/n12. EMD-13921 reads n10 .32 at 45 A.
  - DJ maps get `fold_ft_dj.sbatch out.npz name=map ...` (tree/symfold.py projections, then fold_ft.py next to EMD-13921), reporting |F| n6/n10/n12 per band.
- 20:40 IgM session: SIMPLE symmetry_test can't call 5 vs 6 (EMD-13921 c5 .83 = c6 .83; c2 tops every map, being a flat disc). The pentamer (5 Fc + J chain) has no exact cyclic symmetry; its 5-fold is in the Fab arrangement. DJ maps: the |F| test (fold_ft.py) next to EMD-13921 only.

## 21:00 — |F| test on OUR face-on 2D star classes (fold_ft2d.py, local, seconds)
- Our tpick_all 2D classes 92 / 138 / 246 (the 27k face-on stars), binned to 3.27 A/px:
  - all read n12 .88-.95 at 80 A, n6 ~.47 at 60 A, n12 .51-.59 at 35 A;
  - n10 <= .06 in every band.
- References: EMD-13921 n10 .32 at 45 A; our consensus maps n6/n12.
- By eye at full size (`~/Downloads/igm_ours_faceon_star_classes_92_138_246.png`): 138 is plainly 6-fold (6 holes, 6 outer blobs); 92 and 246 are less regular.
- The earlier 2D ring readout called 92 5-fold, so the two readouts disagree on 92.
- Reading: face-on 2D averages look 6-fold too. That means either face-on hexamers, or 2D alignment putting pentamers in different 60-degree rotations on their pseudo-hexagonal Fc core (the same problem as in 3D). A synthetic pentamer+hexamer set would separate the two.
- 20:25 r300k (54965797): DJ pool scaled to box 200 / 3.27 A by 20:14 (dj_ws2/pool200/dj_pool_sc.simple, all 2,251,109 particles); 300,000 drawn and pruned by 20:23; workers on csbdevel; 6-state solve3D running. Watcher bc93zs8k5 runs fold_ft_dj.sbatch on its 6 state maps at the end.

## 20:45 — user: "both. please proceed autonomously overnight" — two fast tests launched (folder igm_cryolo/ph_tests/)
- Face-on 2D on OUR data, job 54968041 (`face2d.sbatch`):
  - the 27,039 particles of star classes 92/138/246 from pool200 (box 200), pruned (27,039 active at 20:35);
  - solve2D ncls=20 mskdiam=560 in shared memory (32 threads, no worker re-queueing);
  - watcher bti4nxj29 fetches the cavgs and runs the |F| test per class.
- Synthetic pentamer vs hexamer, job 54968346 (`synth.sbatch` + `synth_ph.py`, tested locally first):
  - P = EMD-13921 inside r 140 A, disc normal on z, centred on its pseudo-6-fold axis. J slot = least wedge mass (996 vs Fc 1406-2553).
  - H = C6 sum of a typical Fc wedge not next to J (soft 60-degree partition of unity), scaled to six mean Fc units (mass H/P 1.08).
  - Projections: P is a star with an asymmetric gap, |F| n10 .26 at 45 A (cf. EMD-13921 .32); H is a clean 6-fold star, n6 .63 at 60 A, n10 0.
  - Images: simulate_particles with CTF, then the IgM fixture's quality + colored noise (measure_cn.json). 3D sets use the real tb100k poses (models rotated onto the consensus normal); the face set is face-on +-12 deg.
  - Runs: 3dP 30k pentamers, solve3D 1 state (register-average test); 3dmix 45k P:H 2:1, solve3D 2 states; face 9k P:H 2:1, solve2D 10 classes. Then purity vs truth and |F| per map / class.
  - Watcher bs4n8omyi.
- Still running: DJ class-based 54965177 (watcher b5v7m1hci), DJ r300k 54965797 (watcher bc93zs8k5, runs fold_ft on its maps), IgM session's rep_A/rep_B.

## 21:00 — FACE-ON 2D RESULT (job 54968041, 24 min, 20 classes of the 27,039 face-on star particles)
- |F| test (fold_ft2d.py, the audit/fold_ft.py method on class averages) splits the face-on particles:
  - Pentamers (5 Fab lobes and a visible J-chain gap; |F| n10 .25-.33 at 45 A, cf. EMD-13921 .32): classes 1/7/13/20, 9,441 particles, 8.5-9.3 A.
  - Weaker gap (n10 .14): classes 4/9, 4,373.
  - Regular six-lobed stars (n10 ~0, n12 .7-.8 at 80 A): classes 11/14/18/19, 7,958, 8.5-12 A, as sharp as the pentamers.
  - Junk / tilted / side views: 5,267.
- So 2D CAN resolve the pentamer register for face-on views. The open question is whether the 6-fold classes are hexamers or pentamers whose gap was not resolved.
- Follow-up job 54969247 (`face2d_sub.sbatch`): re-classify the 6-fold set (7,958) and, as a control, the pentamer set (9,441), 10 classes each.
  - If gapped pentamer classes come out of the 6-fold set, those were mis-registered pentamers; if every sub-class stays 6-fold, hexamers are likely.
  - Watcher bu295mvd8.
- Mac ~/Downloads: igm_ours_faceon_reclass20_cavgs_iter029.mrc, igm_ours_faceon_reclass20_pentamer_vs_6fold.png (n10 classes marked D).
- 21:10 first sub-classification (54969247) stopped at its count gate: face2d's solve2D used mkdir=yes, so the classes are in face2d/1_solve2D/face.simple, not face2d/face.simple. Stale copy moved to face_six.simple.unsolved_copy; resubmitted from the solved copy as 54969353 (face2d_sub2.sbatch), watcher bdudzafdh.

## 21:30 — FACE-ON VERDICT: the six-lobed face-on classes are pentamers with unresolved register, not hexamers (calibrated on the synthetic set)
- Sub-classification (54969353, 18 min, 10 classes each):
  - 6-fold set (7,958): stays six-lobed. Sub-classes 2/3/6/7/8/10, 6,897 particles, 9.2-12.8 A: |F| n10 <= .10 at 45 A, n12 .65-.86 at 80 A. Only class 1 (251, 27 A) looks pentameric.
  - Pentamer control (9,441): stays pentameric. Classes 1/2/5/7 (5.3k) n10 .26-.40, gap visible; 2 and 5 show the Fab-1 rod at the gap, as in EMD-13921.
- Synthetic face-on 2D (54968346; 9k known P:H 2:1, 10 classes):
  - purity .914 (chance .662);
  - pure P classes 4/5 (4,408, n10 .19-.26 at 45 A);
  - pure H classes 8/10 (2,117, n6 .77-.80 at 80 A, n10 0);
  - but mixed six-lobed-looking classes 6/7/9 (1,652) are 68-82 % PENTAMERS: 2D puts ~1/5 of known pentamers into register-averaged classes.
- Calibrated test (c6_compare.py: C6-symmetrise each average, 12 A low-pass, r < 170 A, best in-plane rotation with or without mirror):

  | | P vs P | H/S vs H/S | P vs H/S |
  |---|---|---|---|
  | synthetic (H = hexamer) | .991 | .980 | .902 (a real hexamer drops to ~.90) |
  | real (S = six-lobed) | .986 | .979 | .983 (no drop) |

- Verdict: once symmetrised, the real six-lobed face-on classes are the same structure as the pentamers. They are pentamers averaged over 60-degree registers, as the synthetic set shows 2D can produce. No hexamer signal among our ~27k face-on views.
- Caveat: the synthetic hexamer is a construct (C6 of a pentamer Fc wedge, core only, no Fab arms); a real hexamer's difference from C6(pentamer) is assumed comparable.
- Mac ~/Downloads: igm_ours_faceon_subclass_6fold_set_top_pentamer_set_bottom.png plus the two sub-class MRCs.

## 21:40 — user: cls_expansion on one pentamer class
- Job 54971163 (`clsexp.sbatch`, dir ph_tests/clsexp_five5): face-on pentamer-set sub-class 5 (1,837 particles; J-chain gap and Fab-1 rod; n10 .35). Pruned to its members first, avoiding the class= whole-project trap (cls-expansion-class-path-defect). cls_expansion class=5 ncls=20 mskdiam=560 nthr=32 (shared memory). Watcher b9erzf84f fetches cls_expansion_cavgs.mrc.
- 21:30 cls_expansion done (10 s compute): 20 kernel sub-averages of pentamer class 5.
  - pop 61-122, neff 453-805, FRC0.143 8.96-9.62 A (sub 20 12.3 A), repro .32-.50.
  - All are pentamers: |F| n10 .21-.42 at 45 A, with the J-chain gap and the Fab-1 rod.
  - Per-pixel SD across the 20 (each normalised in the 560 A mask): highest at r 110-160 A (.143) on the five Fab lobes and the rod; Fc core r < 60 A quiet (.076). The within-class variability is Fab-arm and rod position, not the core.
  - ~/Downloads/igm_ours_faceon_pentamer_cls5_clsexpansion_20subcls.mrc (user asked for the MRC).

## 22:30 — register fix (user approved: synthetic first, then the 683k consensus); scripts in ph_tests/ (scratchpad reg/)
- register_fix.py: per particle, score the 6 rotations (60k deg) of a pentamer reference about the core's pseudo-6-fold axis at the particle's consensus pose (Fourier slices, CTF, whitened), keep the best, rewrite the pose (R S_j; shift + (R S_j u_j)_xy for an axis through the core centre).
  - Conventions verified against SIMPLE's simulate_particles: cc .998. Plane coefficient (h,k) = Vhat(R^T[h,k,0]), R = euler2m; image = projection shifted by -(x,y); CTF = +(sqrt(1-a^2) sin chi + a cos chi).
- Synthetic test (54971882; 3,000 pentamers, random start registers, fixture noise): accuracy .997 (exact poses), .992 (3 deg pose noise), .994 (20 A reference + 3 deg); 1.000 for the top 75 % by margin; corrected poses = truth + the added noise.
- Docking the reference into the consensus frame (needed for real data):
  - the fixture's maps.json axis is off by 2.6-3.2 deg / 3.4-4.5 px for fine features;
  - band-passed (30-8 A) 60-deg self-correlation of the consensus core is .26 at the refined axis vs .03 at maps.json's;
  - with an FFT shift search the fit peaks sharply (.25-.27 vs ~.09 background);
  - mirror hand favoured at 30-8 A (.267 vs .247) and 12-6.5 A (.112 vs .094): our consensus probably has the opposite hand to EMD-13921.
  - Final docking (dock3.py) puts the reference's own pseudo-6-fold axis on the consensus axis and searches only rotation, hand and height.
- Next: regfix_tb.sbatch on tb100k (100k consensus particles at box 200): score, rewrite poses, reconstruct original vs corrected vs top-half-margin, then the |F| fold check.
- 22:20 coordination: the IgM session (local_126e2179) is implementing the register fix itself (user-approved there: registered reference jreg/2_refine3D_auto in the EMD-13921 frame, box 200; C6-start control jreg/3_refine3D_auto, job 54972007). It asked me not to build it in parallel.
  - Handed over everything in ph_tests/regfix_handover/: scorer, convention check, docking scripts, dockbp/dock2/refdock outputs, and the unsubmitted regfix_tb.sbatch.
  - Told it the mirror-hand finding. dock3 (axis-consistent docking) is still running locally; its output goes to the same folder.
  - My tb100k/683k real-data runs are NOT submitted.

## 22:30 — CORRECTION: the face-on "no hexamers" calibration was too optimistic
- dock3 shows the synthetic P model's true pseudo-6-fold axis is 5.4 px / 2.4 deg from the construction axis (60-deg self-corr .586 vs .033 there).
  - The v1 synthetic hexamer (C6 sum of a wedge about the construction axis) therefore had an artificially different Fc ring.
  - The synthetic register test is unaffected (same axis used throughout).
- Noise-free face-on test (12 A, r < 170 A): C6-averaged pentamer vs a hexamer built about the TRUE axis = .965 (the v1 pair .867).
  - Real six-lobed vs C6(pentamer) .983 (within-set .979-.986) still leans to register-averaged pentamers, but a hexamer sits only ~.02 lower: suggestive, not decisive.
  - The v1 synthetic face-on purity (.91) and the 3dmix run overstate separability; 3dP is unaffected.
- dock3 (axis-consistent docking): mirror hand (12-6.5 A .114 vs .087), phi 14, height 0; slot masses Fc 70-91, J 9.
  - refdock_ref (free fit) is within ~1 px of it (corr .85); my "5 px off" warning to the coordinator was wrong and is corrected.
  - dock3_ref.mrc + dock3.json are in ph_tests/; the coordinator runs the register fix on tb100k with them.
- v2 synthetic face-on test with correct models: P2 = docked pentamer with its axis on z; H2 = C6 of a typical Fc wedge about the true axis, mass 1.15 x P2. Job 54974178, watcher b4qncbv4s.

## 22:45 — RETRACTED: face-on "no hexamers". 2D cannot decide it (synthetic v2, job 54974178)
- Hexamer built about the true axis (ph_tests/synth_face2/{P2,H2}.mrc). Same 9k face-on views, P:H 2:1, fixture noise, solve2D 10 classes.
- Purity .744 (chance .662; v1 with the faulty hexamer .91).
  - Only classes 7/8 (2.3k) are ~90 % pentamer, with a visible gap and n10 .25-.38.
  - The six-lobed-looking classes are P/H mixtures, mostly majority pentamer (56-71 % P); the H-rich class (80 % H) looks identical.
- The real six-lobed face-on classes (7.9k) are therefore undetermined: unregistered pentamers and/or hexamers. Coordinator told.
- Decisive test proposed (synthetic only): 7-way per-particle scoring, 6 pentamer registers + the hexamer, on the same J-chain/gap signal that gives 99 % register accuracy.
  - Models in the consensus frame: ph_tests/synth_ph7/{P2_cons,H2_cons}.mrc (round-trip .998; H C6 self-corr .995).
  - Job 54975432 (16k synthetic, real tb100k poses, P:H 2:1, 3 deg pose noise), watcher bec865o5m.
- 22:55 coordinator approved the 7-way synthetic validation (its real tb100k register run is regfix/tb, job 54974704). Its three asks and my answers:
  1. H2 from the full subunit wedge: already so. The wedge is cut only by angle inside the same 140 A sphere as P2 (Fc + Cmu2 + hinge stub); only the pentamer-specific J chain and Fab-1 rod are absent. Neither model has the flexible Fab arms (EMD-13921 lacks them).
  2. Breakdown by view (face <30, tilted 30-60, side >60), by margin, and the D distribution: new modes synth7d (chunked, per-particle table) + report7 in register_fix7d.py. Smoke test on 40 clean images: all calls right, D P -25..-14 vs H +5..+18.
  3. Per-candidate normalisation kept. The false-H rate comes from the truth-P particles, scored independently of the mix.
- synth7d.sbatch is staged in ph_tests/synth_ph7/; I submit it once job 54975432 has written ph7.mrcs.

## 23:00 — real tb100k register fix made the map WORSE (coordinator's run regfix/tb, job 54974704); my diagnosis
- Coordinator's numbers: FSC 7.43 A (original poses) -> 12.34 A (corrected) / 12.82 A (top half); |F| 45 A n10 .15 -> .01; j counts 20.7k / 28.9k / 20.0k / 4.6k / 6.5k / 19.3k.
- Frames match: core.mrc vs rec_orig corr .947, every flip/transpose <= .12.
- 60-deg self-correlation about the dock3 axis: rec_orig .44 (30-8 A) / .25 (12-7 A); rec_corr .82 / .80. The fix produced a register average.
  - rec_orig is less symmetric than the registered pentamer reference (.59): the tb100k consensus is largely REGISTERED, not a 6-way average. Its 6/12-fold |F| likely comes from the pseudo-hexagonal core.
- The dock3 reference barely prefers any register of rec_orig (j = 0..5: .38 / .42 / .37 / .31 / .31 / .34).
- Hand: against rec_orig the AS-IS hand fits better (.483 vs .421 at 30-8 A; .310 vs .276 at 12-7 A), the opposite of the core.mrc call. A wrong-hand reference makes the calls meaningless. Coordinator told; it should re-dock against rec_orig and score against rotations of the consensus itself.
- 7-way synthetic (job 54975432; correct models, real tb100k views, 3 deg pose noise): P called P .973 (false-H .027), H recall .820, register accuracy .953, overall .922 (chance .664). Detailed breakdown job 54977080. Synthetic success assumes a matched reference.
- 23:05 7-way detailed breakdown (job 54977080, table ph_tests/synth_ph7/tab7.txt):

  | view | false-H | H recall | register acc (P called P) |
  |---|---|---|---|
  | all | .028 | .824 | .953 |
  | face-on <30 (1,351) | .025 | .804 | .953 |
  | tilted 30-60 (5,381) | .030 | .800 | .946 |
  | side >60 (9,268) | .027 | .840 | .957 |

  - Side views are as informative as face-on.
  - |D| > .33 (q25): 75 % called, accuracy .982, false-H .002, H recall .930. |D| > .63: 50 % called, accuracy .997.
  - Register accuracy by margin quartile: .837 / .982 / .996 / .997.
  - Sent to the coordinator with the real-data caveat (matched reference first).

## 23:10 — synthetic de novo 3D (job 54968346, 2 h 12 min; full mask, lpstop 10, box 200, real tb100k views, fixture noise)
- 3dP (30k pentamers, 1 state): the map stays pentameric.
  - Ring n5/n6 at 110 A .68/.02 (max-n5n6 centring; model .70/.00).
  - |F| n10 at 45 A .12 vs the model's .29: partial register mixing, not a 6-fold average.
  - The real tb100k rec_orig reads n10 .15, like this pure-pentamer de novo map.
- 3dmix (45k P:H 2:1, v1 hexamer, 2 states): NO separation, state 1 = 44,997 (purity .671 = chance), state 2 = 3.
  - De novo multi-state failed even with the over-different v1 hexamer, so a missing 6-fold class in the 6-state reps/r300k is weak evidence against hexamers.
  - The per-particle 7-way scoring with a matched reference separates them (.92 vs .66). Coordinator told.
- 23:10 coordinator asked whether P2_cons/H2_cons are matched candidates. Answer: NO (EMD-13921 source, r < 140 A core without Fabs, mirror hand).
  - Built a matched hexamer from OUR consensus: ph_tests/hcons/H_orig_cons.mrc (+ .json). C6 of one subunit wedge of regfix/tb/rec_orig out to r 230 A, about the axis refined on rec_orig (1.15 deg / 0 px from dock3; self-corr .443 vs .441).
  - H_orig_cons C6 self-corr .990; corr with rec_orig .527 at 30-8 A (EMD reference best .42-.48).
  - Caveat: rec_orig wedge masses (r 30-110 A) are 43.9 / 32.5 / 29.0 / 23.9 / 35.8 / 63.9 k, a smooth one-sided pattern, not 5 Fc + empty J (EMD 72/75/9/91/74/70). Either a lopsided map or an axis point ~1-2 px off.
  - Coordinator's tb3 run (job 54977294) scores against rotations of the consensus itself, cross-validated by eo.

## 23:15 — the consensus is REGISTERED (coordinator's cross-validated tb3, job 54977294)
- Scoring against the consensus's own six rotations (eo=0 vs the odd half, eo=1 vs the even half): j counts 97,549 / 937 / 317 / 160 / 204 / 833.
  - Confident half: 49,964 j=0 vs 38 other.
  - FSC 7.43 A (original) / 7.60 A (all calls) / 7.43 A (confident); 60-deg self-corr unchanged.
  - So the tb100k/683k consensus is >= 97.5 % registered and needs no register fix. Its 6/12-fold |F| comes from the pseudo-hexagonal core, not register averaging (consistent with the synthetic 3dP: a pure-pentamer de novo map has |F| n10 .12 vs rec_orig .15).
- Next on the coordinator's side: a real-data 7-way P/H call, with the consensus half-map rotations plus a matched hexamer.
  - Built ph_tests/hcons/H_orig_dock3axis.mrc: C6 of a full subunit wedge (r 230 A) of rec_orig about the exact dock3 axis, template slot 0, J-like slot 3 unused; C6 .990, corr rec_orig .524, .980 with H_orig_cons.
  - Caveat sent: uneven slot masses 1.13 / .75 / .65 / .52 / .99 / 1.96 x mean; equalising would need a ~20 A shift, contradicted by the centred hub. The axis is the least-validated piece; tb3's j=0 dominance can't validate it.
- 23:25 axis-free register check (coordinator's ask) running as job 54977934 (ph_tests/facecons/face_cons_avg.py). It averages the face2d six-lobed (11/14/18/19) and gapped (1/7/13/20) particles with their 683k CONSENSUS poses (t2more_it3, via pool rows from face2d/states.txt): pose-matrix in-plane mapping (handles top/bottom mirroring), shift /3, phase flip. Control: the synthetic face-on P/H set averaged with true poses. Transform checked locally: synthetic top and bottom views map onto the top-view projection equally (corr .844-.849). Watcher bqer54fby.

## 23:35 — real-data pentamer/hexamer call on tb100k (coordinator asked me to run it; job 54978042, ph_tests/real7/)
- Cross-validated: eo=0 rows (regfix/tb3/rows_e0.tsv) vs the six dock3-axis rotations of the ODD half map plus H_odd; eo=1 vs the EVEN half plus H_even.
  - P_* = rec_orig half maps under a 230 A soft sphere.
  - H_* = C6 of each half's slot-0 subunit wedge about the dock3 axis: C6 .990, corr .994 with the full-map H.
  - The uneven slot masses are reproduced in both halves (1.14 / .85-.86 / .76-.78 / .64-.65 / .93-.95 / 1.65-1.67), so they are real structure, not noise.
- Then: D distribution vs the synthetic calibration; counts called H at D > 0 / .33 / .63.
- Subsets N = max(2000, #called H): most H-like vs strongest P, via state_patch.py on copies of tree/tb100k/igm_sub_sc.simple (poses untouched, row-for-row pose assert).
- reconstruct3D (c1, mskdiam 560, mkdir=no as the coordinator did), then slot masses / J:Fc / C6 / |F| for full and half maps. Watcher in place.
- 23:40 face-on consensus-pose averaging: synthetic control OK (truth P: 60-deg self-corr .52, |F| n10 .35 at 45 A; truth H: .92, n10 0; top and bottom views agree). The real part failed on too many open files (17,905 stacks memmapped); reader fixed to cache headers only; rerun job 54978096 (face_cons_avg2.py) on the existing dumps, watcher bdq1tmyke.

## 23:45 — REAL tb100k PENTAMER/HEXAMER CALL: no hexamer population (job 54978042, 4 min scoring)
- D = hexamer minus best register (cross-validated half maps), quantiles 5/25/50/75/95:
  - real -3.31 / -2.15 / -1.49 / -.96 / -.35: all on the P side, no H mode.
  - synthetic P -1.82 / -1.20 / -.83 / -.50 / -.09; synthetic H -.27 / .10 / .32 / .56 / .92.
- Called H: D > 0 1.2 %, D > .33 .22 %, D > .63 .04 %; equal to or below the synthetic pentamer false-H rates (2.8 % / .2 %).
  - By view: face-on 4.0 % (synthetic false 2.5 %), tilted .8 %, side 1.1 %.
- Registers 97.8 % j=0 (confirms tb3: the consensus is registered).
- Decisive: the 2,000 most H-like and the 2,000 strongest-P particles reconstruct alike (FSC .143 ~26-27 A):
  - subH slot masses 1.13 / .80 / .84 / .73 / .86 / 1.64, J/Fc .70, C6 .565
  - subP 1.11 / .87 / .76 / .72 / .91 / 1.62, J/Fc .69, C6 .569
  - half maps agree; a hexamer subset would give equal slots and C6 ~.9.
- Reading: no detectable hexamers inside the consensus set (< ~1 %). Caveat: the test only covers particles in the consensus.

## 23:55 — face-on particles averaged with their CONSENSUS poses (job 54978096): the gap is lost
- All three groups give sharp, identical-looking averages (ring aligned; top and bottom views agree), so poses and shifts are accurate:
  - gapped 1/7/13/20 (9,235): 60-deg self-corr .21, |F| n10 .04 at 45 A
  - six-lobed 11/14/18/19 (7,667): .25, n10 .00
  - weak-gap 4/9 (4,278): .23
  - (self-corr is low partly because rotation was about the image centre, not the axis.)
- Their 2D class averages had n10 .25-.40, and the synthetic pentamers averaged with true poses give n10 .35. Aligned by consensus poses the 10-fold signal vanishes and the |F| is 6/12-fold (n12 .68-.82 at 80 A).
- Suggests the consensus does not keep the face-on particles' gaps in register, even though the bulk (mostly side/tilted views) scores as registered.
- Decisive follow-up: per 2D class, the histogram of consensus in-plane angle minus 2D angle (one peak = registered; 60-deg comb = scrambled). Job 54978217.

## 00:05 (Nov 6) — face-on views are partly misregistered in the consensus; side/tilted registered
- Consensus in-plane angle minus 2D angle per GAPPED class (the gap pins the 2D angle; six-lobed classes are ambiguous by 60 deg and are skipped):
  - 13 top (2,625) and 20 bottom (2,317): ~90 % in one peak.
  - 1 bottom 46 / 31 %, 7 top 47 / 42 %: two peaks 60 deg apart. The small subsets are three-way splits.
  - Weighted: ~68 % in the main register, ~32 % off (mostly +-60).
- real7 register calls by view: j != 0 at <15 deg 22.8 % (n 1,930; mostly j = 1/5); 15-30 deg 8.8 %; tilted 1.5 %; side 1.3 %.
  - Face-on margins are low (median .76 vs 1.4-1.5). Confident half: <15 deg 5.9 %, others ~0.
- Reading: the side/tilted majority (~92 %) is registered and sets the map. ~20-30 % of near-face-on particles sit one register off, which is why the face-on 2D classes look six-lobed (with the synthetic v2 result). Too few to move the FSC.
- No change to the hexamer answer: none detectable inside the consensus. Coordinator told.

## 00:20 — coordinator: "registered" and "no hexamers" are NOT established
- Reason: the user sees a hexamer in the 3.85 A map; our consensus J/Fc is .59-.70 vs ~.11 in EMD-13921; tb3/real7 scored against consensus-derived references (self-referential).
- Two asks, both running:
  1. Reference-free per-particle registers for the 9.4k gapped face-on particles from their 2D gap + consensus pose (face_register.py, job 54978407).
     - Step A: within class/side, bring each particle to the majority Q.
     - Step B: align class gap angles in the map frame.
     - Validation: the pooled consensus-pose average must regain the gap. Output: corrected poses (box 600) for patching.
  2. 3D fix with an EXTERNAL clean-gap reference (job 54978511, ph_tests/asis3d/).
     - Reference: EMD-13921 in its own hand (as-is), docked on the consensus axis against rec_orig. Fit .501 (mirror was .421).
     - 60k-rotation fits .501 / .454 / .370 / .369 / .376 / .424: one dominant register plus +-60 pull.
     - Reference J slot at the consensus's lightest slot (~196 deg), J/Fc .12; mask = +-45 deg wedge about the J slot, r 15-170 A.
     - Score tb100k with the masked reference, rewrite poses, reconstruct all / confident / low-margin (bias control). Judge by J/Fc in half maps vs the original .59.
- 00:25 coordinator: "side/tilted registered" and "no hexamers" stay UNPROVEN.
  - If ~92 % were J-registered, the consensus J slot would be near EMD's, not .59-.70; the face-on ~30 % of ~8 % adds only ~2-3 %.
  - The side/tilted j != 0 rates (1.3-1.5 %) came from consensus-derived references, so they show self-consistency, not registration.
  - My two-population estimate: J slot = .12 r + 1.0 (1 - r) = .59 gives r ~ .5, so only about half registered.
  - Coordinator's decisive runs use a reference with NO gap (EMD-13921 as-is, J sector filled by the neighbouring Fc rotated 60, Fab-1 rod kept, 30 A low-pass; jreg/emd_jfill_b200.mrc): jreg_jfill (54978458, 150k) and fullreg (54978459, 1.64M). An empty J slot in both half maps would make the gap data.
  - Its metric audit/jslot.py (EMD frame): EMD .24, first jreg .42, C6 1.0.
  - Complementary to my asis3d (external clean-gap reference + low-margin bias control) and facereg (face-on, reference-free).

## 00:35 — FACE-ON REGISTERS (reference-free; job 54978407): the gap is restored
- 9,235 gapped face-on particles. Registers relative to the consensus pose, j = 0..5: 6,265 / 934 / 108 / 216 / 337 / 1,375; 32 % change, mostly +-60.
- In-peak before correction: 13 top and 20 bottom .91; the other subsets .35-.55.
- Pooled consensus-pose average, before -> after:
  - 2D J/Fc .22 -> .11
  - slots/mean 1.20 / 1.21 / .86 / .26 / 1.09 / 1.39 -> 1.19 / 1.17 / 1.09 / .13 / 1.18 / 1.25
  - |F| n10 at 45 A .04 -> .14
- The common gap angle 198.4 deg matches the consensus's lightest slot (197.7) and the as-is EMD J slot (195.8): three independent sources agree.
- Corrected poses (box 600, t2more_it3 rows): ph_tests/facecons/facereg_registers.txt + facereg_before_after.mrcs. Sent to the coordinator.

## 00:45 — 3D fix with the EXTERNAL EMD reference FAILED (asis3d, job 54978511)
- Calls j = 0..5: 13.4k / 14.0k / 31.6k / 21.7k / 9.6k / 9.8k (13 % stay; 120-180 deg dominate).
- FSC 12.3 / 13.6 / 15.2 A (all / confident / low) vs 7.43 original.
- J/Fc 1.10 / 1.11 / 1.09 (original .59); 60-deg self-corr .85 / .79 / .88. The J slot FILLED and the low-margin control looks the same: systematic bias.
- Cause: EMD's non-gap asymmetry, a heavy slot at J-60 (Fab-1 rod root, 1.39 x mean), vs ours light at J-60 (.76) and heavy at J+120 (slot 5, 1.67). The masked score matched mass, not the gap. Coordinator warned that its jfill reference keeps the rod too.
- Next: GAP-ONLY reference = our consensus C6 Fc ring (H_orig_dock3axis) minus the J-slot wedge, cylinder r < 125 A; slots/mean 1.20 / 1.20 / 1.11 / .21 / 1.07 / 1.20, J/Fc .18.
  - Validate on the 9,235 face-on particles with reference-free registers, then tb100k (all / confident / low; half-map J/Fc).
  - Job 54978719 (ph_tests/gaponly/).

## 23:55 (2026-11-05, Mac clock) — GAP-ONLY reference (job 54978719): passes the face-on check, but the tb100k fix is only partial
- Clock correction: the headers above from the "23:55" face-on consensus-pose entry through "00:45" are mislabeled. Scheduler end times on 2026-11-05: face_cons (54978096) 23:21, facereg (54978407) 23:37, asis3d (54978511) 23:40, gap-only (this entry) 23:50.
- Reference: ph_tests/gaponly/G_gaponly.mrc + gaponly.json (box 200, 3.27 A, consensus frame, dock3 axis). Slots/mean 1.20 / 1.20 / 1.11 / .21 / 1.07 / 1.20, J/Fc .18.
- Face-on validation against the 9,235 reference-free registers:
  - Agreement .783 (chance ~.17; "all j = 0" .678).
  - By margin: >= q25 .864, >= q50 .908, >= q75 .932.
  - Among the 2,970 particles the face-on step moved: .652.
- tb100k calls j = 0..5: 32.3k / 23.8k / 2.7k / 3.1k / 17.4k / 20.7k (68 % change). Margin quartiles .24 / .54 / .97.

| Map | FSC .143 | J/Fc (half maps) | Fc CV | 60-deg self-corr |
|---|---|---|---|---|
| original | 7.43 A | .59 / .60 | .31 | .44 |
| all | 9.91 A | .46 / .48 | .08 | .82 |
| confident half | 13.08 A | .43 / .43 | .08 | .78 |
| low-margin half (control) | 13.08 A | .51 / .55 | .09 | .85 |

- Reading:
  - The J slot empties only partway, not to the .1-.2 target.
  - Confident and low-margin barely differ, so part of the shift is the reference's own C6 ring pulling the map (Fc slots equalise, self-corr rises, resolution drops).
  - The face-on calls are reliable; the side/tilted calls are not.
- Sent to the coordinator with two cautions:
  1. A gapped reference can't serve as the decisive refinement reference (a gap in the output would be ambiguous). Use it only as a per-particle scorer on face-on/confident particles.
  2. Its jfill reference keeps EMD's rod (heavy at J-60; ours heavy at J+120). Rod-to-heavy-sector matching would put a data-driven gap opposite EMD's J (EMD J + 180). Check that sector in the six-sector readout.

## 00:20 (2026-11-06) — six-lobed face-on classes are mostly GAPPED PENTAMERS (job 54979324; coordinator's design)
- Set-up:
  - 7,667 six-lobed (11/14/18/19) and 9,235 gapped (1/7/13/20) face-on particles, scored with the gap-only reference G at their consensus poses.
  - Per-particle controls at each particle's own pose, CTF and signal level, with background noise: C6 ring (no gap) and G at a random register, each at the consensus pose and perturbed by 3 deg / 1 px.
  - Cross-fitted averages: the register called from one Fourier half (kx+ky parity) places the other half.
- Margins. Median, and fraction >= the gapped median:
  - gapped 1.02 / .50
  - six-lobed .85 / .41
  - C6 control .61 / .27
  - pentamer control 1.5 / .73
  - Signal levels match.
- Gap excess vs the gapped set: .66-.82. Even/odd half-call agreement .42 vs .51 (C6 .17-.20), which gives .73. Every six-lobed class/side .48-.60 (gapped .52-.73).
- Cross-fitted J/Fc at the J slot (no gap = 1):
  - all: six-lobed .40, gapped .21
  - margin >= .8: .27 vs .11
  - margin >= 1.3: .24 vs .09
  - The gap sits at the J slot in all of them.
- Reading:
  - >= ~75 % of the six-lobed particles carry a gap as visible as the gapped set's.
  - <= ~25 % (~1,900, ~11 % of the 16.9k) show none: hexamers or gap-invisible pentamers, which this test can't separate.
  - Six-lobed 2D classes are not hexamer evidence.
- Why the cross-fit is clean: the final placement (best of six rotations) is the same whatever register the consensus chose. The consensus-pose six-lobed average (J/Fc .43) is not used as evidence, because consensus registers can be noise-fitted.
- CORRECTION: 2D slot windows in face_register.py took each average's own six-fold phase with a sign error, so they drifted between averages. With fixed windows at the 3D slots (17.75 + 60k; J 197.75), the gapped face-on before -> after is .16 -> .08 (reported .22 -> .11).
  - The 3D stored windows are right: the 2D gap removes exactly the 180/210 sub-lobe pair. EMD as-is .13, consensus .60 unchanged.
- Files: ph_tests/sixlobe/ (sixlobe_avgs.mrcs + .txt + .png, sc_*.txt, rows.tsv); logs igm_sixlobe.54979324.out, igm_sixlobe_rep.54979384.out. Sent to the coordinator.
- csbdevel is at the account cap (632/642 running, 37 pending), so this ran on norm (the pending job was moved; see the partition memory).

## 01:35 (2026-11-06) — DJ class-based solve3D states: |F| test and J gap (job 54980662, `dj_ws2/checks_3D_dj/`)
- The solve3D step of job 54965177 split DJ's 409,710 core side-view particles into three states:
  - s1 187,017; s2 87,937; s3 134,756.
  - FSC 0.143 is 16.35 A for each state.
  - The tests ran on 3x-binned copies (3.27 A/px).
- |F| test (centring-free), largest n10 share per state:
  - s1 .29 at 45 A, .31 at 28 A;
  - s2 .26 at 45 A, .19 at 28 A;
  - s3 .45 at 35 A, .26 at 80 A.
  - References: EMD-13921 .32 at 45 A. Our consensus maps reach at most .16 in any band; n12/n6 dominate them.
- J gap (`regfix/jslot_self.py`, lightest sector over the mean of the other five):
  - s1 .36;
  - s3 .42;
  - s2 .06, unreliable: the sector next to it is also light (.26), the centre sits at the +-8 px search edge, and ring self-corr is .32.
  - References: EMD-13921 .13 in the same run; our consensus .60-.63.
- Disc-normal projections (`checks_3D_dj/dj3d_proj.png`):
  - all three are blurred by rotation about the normal and show radial streaks (the core set is side views only);
  - s3 shows a ring of beads;
  - s2 shows a core with one rod, like EMD-13921's Fab-1 rod.
- Reading: all three states lean pentamer by both tests, unlike our consensus. None reads 6-fold. At 16 A with rotational blur this is suggestive, not decisive. refine3D (the brief's last step) would sharpen the J-gap call.
- Maps: `~/Downloads/igm_dj_classbased_solve3D_state{1,2,3}_{187k,88k,135k}_16A_3.27apix.mrc`.
- The addon on the 701k super set is in its last stage. Watcher bgmhma2i9 runs `checks_3D_dj/statechk.sbatch` on the addon maps when it ends.

## 03:30 (2026-11-06) — CORRECTION: the DJ class-based states are at the C6 level (job 54982827, `checks_3D_dj/regemd_dj.sbatch`)
- Prompted by the coordinator's 03:15 readout trap: `regfix/jslot_self.py` takes the disc normal from the whole thresholded map, so the arms set it. `tree/symfold.py`, which the |F| n10 test projects through, does the same.
- The 01:35 "all three lean pentamer" call (J gap .36/.42, |F| n10 .26-.45) is withdrawn.
- Deposited-frame test (reps/regemd protocol):
  - each 3x-binned state map and its mirror was docked into EMD-13921 (dock_volpair lp=20 mskdiam=280 trs=15);
  - `regfix/regemd.py` read them on EMD's core axis and centre, alongside the coordinator's control dockings.
- Register contrast (core CC at k=0 minus the best CC at k!=0):
  - controls: EMD +.189, EMD J-filled +.075, EMD-C6 +.001, our consensus -.020;
  - s1 +.033 (mirror +.022);
  - s2 +.034 (mirror +.012);
  - s3 -.042 (mirror -.003);
  - all within the 12 replicate classes' -.034..+.038, i.e. at the C6 level.
- Sectors read J-first are gradients or irregular, not EMD's single empty sector (s1 .19 .28 .69 1.15 1.83 1.04).
- Core CC vs EMD at k=0: .67 / .63 / .77 (mirrors .66 / .57 / .75), below EMD-C6's .895. As for the coordinator's maps, "register-averaged" is not separated from "different near the core". Here the side-view-only rotational blur about the normal also erases register contrast.
- The maps in ~/Downloads stay as maps; they are not pentamer evidence.
- Next:
  - the addon watcher (b12ickf0t) now runs this regemd test (not statechk) on the addon maps when the addon ends;
  - refine3D on the addon states is held for the coordinator or the user to decide;
  - when r300k (54965797) ends, its maps (box 200) get the same regemd test via symlinks named r300k_state0N_b200.mrc.

## 04:00 (2026-11-06) — addon done (54965177 COMPLETED, 8 h 30 min); addon maps also at the C6 level (job 54982994)
- solve3D_addon of the 99-class super set (701,487 particles) onto the frozen core states:
  - state 1: 264,029 (core 187,017 + cohort 77,012);
  - state 2: 210,550 (87,937 + 122,613);
  - state 3: 226,908 (134,756 + 92,152);
  - 1,549,622 deselected.
  - The log warned that sampling coverage was short: the least-visited unit got .22 visits per particle; one sweep needs 56 draws, 12 were planned.
- Deposited-frame register contrast (same protocol as 03:30; controls EMD +.189, J-filled +.075, C6 +.001, consensus -.020):
  - state 1 +.034 (mirror +.025);
  - state 2 +.029 (mirror +.009);
  - state 3 -.044 (mirror -.007).
  - Core CC vs EMD-13921 at k=0: .68 / .62 / .77.
  - The addon maps are nearly the core-state maps, and the sectors are unchanged.
- All six DJ class maps (core and addon) are at the C6 level. The class-based branch of the brief ends here.
- refine3D is held until the coordinator or the user decides.
- r300k: watcher bknpxqwtp replaces bc93zs8k5, which would have run the whole-map-normal |F| test. When 54965797 ends, it links the box-200 state maps as `checks_3D_dj/r300k_state0N_b200.mrc` and runs `regemd_dj.sbatch` on them.

## 05:45 (2026-11-06) — r300k done (54965797 COMPLETED, 9 h 53 min); state 6 is the first map above the C6 level (job 54983473)
- r300k = random 300k of DJ's extraction at box 200 / 3.27 A, solve3D 6 states, full mask, no 2D selection.
- Final populations and FSC 0.143 (`dj_ws2/pool200/r300k/solve3D.log`):
  - s1 102,332 (16.0 A), s2 34,432 (16.0), s3 28,129 (16.0), s4 33,548 (15.2), s5 28,901 (14.9), s6 72,658 (16.0).
- Deposited-frame register contrast (same regemd protocol; controls EMD +.189, J-filled +.075, C6 +.001, consensus -.020):
  - s1 +.008 (mirror -.008), s2 +.026 (+.021), s3 +.023 (.000), s4 +.028 (.000), s5 +.011 (+.005) = C6 level;
  - s3 and s5 are junk-like (core CC vs EMD .36 / .18);
  - **s6 +.138 (mirror -.105)**: core CC at k*60 = .740 .571 .539 .602 .597 .490, a SHARP k=0 peak, not the smooth n=1 gradient the C6-level maps show;
  - s6 sectors J-first .31 1.03 1.29 .89 .59 1.20 (EMD .12 1.02 .93 .90 .95 1.20): J sector partly empty, second dip at k=4 (.59) = possibly a register mix;
  - the mirror's best is k=3 (.706) < as-is k=0 (.740), so the as-is hand fits better, as for the consensus.
- s6 sits between EMD J-filled (+.075) and EMD (+.189), outside every other map read so far (-.044..+.038, incl. the coordinator's 12 replicate classes).
- NOT yet verified: solve3D saved no half maps, so noise is not ruled out. The decisive check is a rec3D of s6's particles with their own poses, even/odd, then regemd on each half (both halves > +.075 at k=0 = real). Not launched: it needs a go-ahead.
- Maps at the top of ~/Downloads: `igm_dj_random300k_state6_73k_16A_registerhit_EMD13921frame_3.27apix.mrc` (s6 docked into EMD-13921) and `igm_EMD13921_deposited_pentamer_box200_3.27apix_overlay_reference.mrc` (same frame, for overlay).
- Note: the class-based (2D-selected side-view) runs on DJ's data gave no register; the random subset, which keeps face-on and tilted views, gave one class with it.

## 06:00 (2026-11-06) — r300k state-6 half-map check PASSES the criterion (job 54983521, 1 min 54 s); halves are not independent
- Coordinator OK'd the check. reconstruct3D `nstates=6 state=6` (selection by the state label, no import_particles) on a copy of the final
  r300k project: `dj_ws2/pool200/r300k/2_s6_halves/1_reconstruct3D/`, gridding, box 200 / 3.27 A, postprocess=no.
  - cFAR .8630 = solve3D's final state-6 value: the same particles, poses and even/odd split give the same reconstruction.
- Register contrast | core CC vs EMD at k*60 deg, k=0..5 | sectors J-first:
  - even +.138 | .739 .571 .538 .601 .596 .489 | .32 1.02 1.28 .90 .60 1.20
  - odd  +.137 | .741 .574 .541 .604 .596 .491 | .30 1.00 1.29 .91 .60 1.20
  - full +.138 | .740 .572 .539 .603 .597 .490 | .31 1.03 1.28 .89 .59 1.20
  - solve3D's own halves: even +.138, odd +.137, same profiles.
  - Criterion (both halves > +.075, k=0 best): met.
- solve3D aligned both halves against ONE shared reference in all five stages:
  - the stage manifest has l_lpset=yes in every stage, and filt_mode was nonuniform_lpset with the matching lp promoted;
  - so simple_matcher_refvol_utils (same in the master1030 build that ran) sets l_use_merged_nu_ref and copies the merged recvol_state06_nu_filt.mrc into both half references.
  - The halves therefore agree by construction at 20 A. The check rules out reconstruction noise, not a feature the shared reference imprinted during alignment.
- Correction to 05:45: solve3D did save half maps (recvol_state06_even/odd.mrc); `ls | head` hid them behind the per-iteration FSC PDFs.
- No refine3D run (per coordinator).

## 07:45 (2026-11-06) — s6 as the MATCHED scorer reference on tb100k (coordinator's ask; job 54984845, `ph_tests/s6score/`)
- Coordinator's independent check runs separately: two disjoint random 300k draws of the box-200 pool (rB 54984558, rC 54984559, `reps/djrep/`), same solve3D protocol, regemd on every state.
- This job:
  1. s6 (r300k `rec_final_state06` + solve3D's even/odd halves) pre-docked onto P_asis (EMD as-is on the dock3 axis) with dock_volpair; s6's own pseudo-6-fold axis moved onto dock3; phi (full turn) and height searched against rec_orig for both hands; the as-is hand is used (`dock_s6.py`).
  2. asis3d's J/gap/rod mask (wedge +-45 deg about s6's own J slot, r 15-170 A).
  3. Face-on check vs the 9,235 reference-free registers, each s6 half separately.
  4. tb100k cross-validated as in tb3: eo=0 rows vs the s6 ODD half, eo=1 vs the EVEN half. Gaussian lp 16 A (s6's resolution); scorer band 150-12 A as before.
  5. Rewrite poses on copies of the tb100k project; reconstruct all / confident / low-margin; FSC; J/Fc in both halves in the consensus J window (197.75); j and margins by view.
- Coordinator's success bar: face-on agreement >= .85 (gap-only .783); FSC not worse than ~7.6 A (orig 7.43); J/Fc <= ~.35 in both halves (consensus .59/.60).
- 07:56 first attempt (54984845, `ph_tests/s6score/`) CANCELLED after docking:
  - dock_volpair onto P_asis was unreliable: full .723 and even .722, but odd .878; the odd half sat in another orientation (corr .09 vs full);
  - s6's own axis was 4.7 deg / 9.9 px off and the fit to rec_orig only .135 (EMD as-is .501), with the J slot at 114 deg;
  - its face-on numbers (.31 / .23) are meaningless; outputs kept, not used.
- 08:10 v2 (54985371, `ph_tests/s6score2/`, `dock_s6b.py`): starts from s6 already docked into EMD-13921 by regemd (full + solve3D halves).
  - The own axes of EMD and s6 are put on dock3 in both senses and both hands; phi and height are searched against rec_orig.
  - Gates before scoring: s6 halves vs full >= .90 in the EMD frame; EMD placed the same way reproduces P_asis (>= .90); s6's axis within 5 deg / 4 px of EMD's.
  - Composition helper checked locally: .995 vs sequential application (a 1 px error gives .90).
- 08:13 v2 (54985371) stopped at its gates; nothing was scored.
  - In EMD's frame, s6 vs EMD-13921 correlates only .06 at 30-8 A, while s6's halves agree .99 there. Its dock_volpair cc .877 is a 20 A-or-coarser match.
  - So s6's 30-8 A content is shared between the halves but not EMD-like: shared-reference imprint, in line with its solve3D halves not being independent.
  - It also left s6's own 60-deg axis unfindable at 30-8 A (self-corr .10).
  - The EMD control placed EMD up/down/hand-flipped against rec_orig at best .464, vs P_asis .501; similarity to P_asis only .716. Method issue (full-map own axis); not pursued.
- 08:20 v3 (54986097, `ph_tests/s6score3/`, `dock_s6c.py`):
  - Anchor = v1's odd half docked onto P_asis (cc .878, the same as s6 vs EMD in EMD's frame); the full map and even half are docked onto the anchor.
  - Gates: halves vs full >= .95 at 40-16 A; s6 in EMD's register on dock3 (core CC vs P_asis rotated 60k: k=0 best).
  - The scorer's low-pass comes from s6's shell correlation with rec_orig (first shell < .5, clamped 16-30 A).
  - Hand: full-turn phi/height vs rec_orig at 40-16 A and 60-20 A.
  - Mask about s6's lightest consensus slot window.
- 08:28 HOLD: the coordinator asked to hold every s6 scorer run until its band check (job 54986056, regfix/regband.py) is back.
  - v3 (54986097) was cancelled at 56 s, before any scoring. Its pre-scoring checks had already printed (`igm_s6score3.54986097.out`).
  - Placement works: full and even halves docked onto the odd-half anchor at cc .998; halves vs full .97 at 40-16 A.
  - In the consensus frame, core CC vs P_asis rotated 60k = .740 .578 .542 .606 .597 .501 (contrast +.134) reproduces regemd's +.138.
  - Shell CC, soft core mask; s6 vs rec_orig | vs P_asis (EMD) | even vs odd:

    | band (A) | vs rec_orig | vs P_asis (EMD) | even vs odd |
    |---|---|---|---|
    | 120-60 | .575 | .776 | .999 |
    | 60-40 | .198 | .427 | .996 |
    | 40-30 | .069 | .209 | .990 |
    | 30-25 | .039 | .092 | .954 |
    | 25-20 | .051 | .059 | .906 |
    | 20-16 | -.012 | -.026 | .815 |
    | 16-9 | ~0 | ~0 | .84-.90 |

  - Reading: s6 matches EMD and the consensus only at >= ~40 A, while its halves agree down to 9 A: shared-reference imprint. Its 16 A FSC is not real resolution, and the +.138 register contrast is a >= 40 A feature.
  - The ready placement for a scorer run, if the band rule passes: `s6score3/anc_{full,even,odd}.mrc`. The scorer low-pass would be 30 A, the clamp; the signal says >= 40 A.
- 08:45 s6 RETIRED as the register reference by the coordinator (band check 54986056):
  - 60-30 A: s6 k=0 .26 < C6-level s1 .30; its relative 60-deg drop .62 > EMD's .37.
  - 30-20 A: every de novo map ~0.
  - All s6 scorer and docking work stopped; no new IgM direction until the user picks the route.
  - Nothing of mine runs. The running jobs on the account are the coordinator's djrep replicates (54984558/9), left alone.
  - Correction sent: the "J 158 deg off, fit .149" numbers came from v2's misdocked placement (failed gates), not a valid measurement. The verdict rests on the band check alone.

## 12:20 (2026-11-06) — register work CLOSED (user decision, relayed by the coordinator)
- IgM tree input = the 683k consensus as is (t2more_it3.simple), with the register caveat. RELION symmetry-expansion pilots v1/v2 failed for side/tilted views (see the coordinator's 12:10 entry, symexp/).
- This session stands down: no further IgM jobs or analysis. Nothing of mine is running on the cluster.
