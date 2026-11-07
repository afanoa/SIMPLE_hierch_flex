# IgM (EMPIAR-11077) overnight processing report

Working dir: `/mnt/beegfs/elmlund/afan/igm_cryolo/` on batch.ncifcrf.gov. SIMPLE build `~/SIMPLE_master1030`.

## Stage 1: picking (done 2026-11-02)
- crYOLO 1.9.6, fine-tuned from the general model on 21,257 template picks from 996 micrographs.
  - Picks came from SIMPLE reference picking with 21 templates: 2D classes plus reprojections 32/38/47 of EMD-13921.
  - 2D-cleaned: flat classes 46, 54 and 24 side-view classes. Best validation loss 0.0383 at epoch 5.
  - Weights: `tpick/cryolo_train/igm_ft.h5`.
- Full set at threshold 0.3: **1,659,513 picks** over 43,331 micrographs.
  - Picks per micrograph: median 37 (p5 12, p95 68). 135 micrographs have none.
  - At threshold 0.1: 2,665,622 picks, kept in `full_ft/run_NN/out/CBOX`.
- Box files: `full_ft/box_thr0.3/` (320 px). Pilot overlays: `~/Downloads/igm_cryolo_trained*/`, `~/Downloads/igm_cryolo_full_picks/`.

## Stage 2: crYOLO extract + solve2D: rejected (2026-11-03)
- Extracted 1,498,284 particles at box 600 (job 54896224). solve2D at iterations 16-21 showed rods and blobs but no 5-fold pentamer class. The user agreed.
- A/B on the same 1,000 micrographs (job 54897578): crYOLO picks gave no star class at iteration 11. The job timed out waiting for workers. Template picks with the same pipeline gave 5-fold classes #46 and #54. Cause: crYOLO training was ~95% side views.
- Cancelled master 54896224 and workers 54907554-73 by ID. Outputs kept in `ft_s2d/`.

## Stage 3: full template pick (job 54897549) -> extract + solve2D
- SIMPLE `pick` with the 21 templates in `tpick/pickrefs_v3.mrcs`: 18 real class averages plus reprojections 32/38/47. 80 tasks over all 43,331 micrographs, output in `tpick_all/pick/`.
- Then `tpick_all/all_s2d.sbatch`: extract at box 600, solve2D with ncls=300 and mskdiam=560, 20 workers.
- Gate: no 3D until the full 2D shows pentamer (5-fold) classes.
- Extraction done 08:39 (20 workers, 43,331 stacks; 2,400,639 particles active in 2D). solve2D at stage 4 (iter 31, `~/Downloads/igm_tpickall_cavgs_iter031.*`) shows 3 full 5-fold star classes plus many small 5-fold core classes: gate PASSED. Letting solve2D finish before class selection.

## Stage 3b: in-plane rotated templates (test, 2026-11-03)
- Found the picker never rotates templates, and `pickrefs_v3` skipped `make_pickrefs`, which adds the rotations. Built `tpick/pickrefs_v4.mrcs`: the 21 templates rotated in 15-degree steps, 504 in total.
- Re-picked the same 1,000 micrographs (job 54911883, 11 min per 50 micrographs), then extract + solve2D with ncls=100 (chain 54911882).
- Full-star (flat) classes: v4 classes 63 + 67 = 1,276 of 51,571 (2.5%) vs v3 classes 46 + 54 = 1,166 of 55,306 (2.1%). No material gain: template rotation was not what limited flat views. Classes in `~/Downloads/igm_1000mic_rotrefs_cavgs_iter030.*`.
- Full solve2D finished at iter 31. Full-star classes 92, 138 and 246 hold 27,039 particles. My selection in the picker keeps 190 classes (1,628,358 particles, 67.8%) and drops 110 junk classes (ice discs, carbon edges, gradients, blobs, empty). Awaiting the user's review.

## Planned 3D workflows (cheap tests first)
- T0: solve3D_cavgs on the kept class averages, c1.
- T1: solve3D c1 on a 150k random subset of the selection.
- T2: strict-core solve3D, then solve3D_addon with 100k dropped + 100k kept particles; the addon_diag cohort map tests the dropped classes.
- T3: solve3D nstates=3 on the T1 subset.
- Afterwards: refine3D_auto on the winner with all its particles.
- Launched 2026-11-03 in `tpick_all/3D_tests/`, all derived from `igm_all.simple`. Counts verified: sel190 = 1,628,358, sub150k = 150,000, core15 = 192,481, t2super = 392,481. Jobs: T0 54918472 (solve3D_cavgs), T1 54918474, T2 54918476 (core solve3D, then solve3D_addon), T3 54918478 (nstates=3).
- The first prep (54916945) ran selection without mkdir=no, so the projects were unchanged. Its four masters were cancelled by ID before any workers started; outputs moved to `3D_tests/bad_prep/`.
- T0 finished: its map from class averages is poor, by the user's judgement (13.3 Å). T4 (54921687) added: refine3D_auto on the same 150k, with EMD-13921 low-passed to 30 Å as vol1 and ref_pose_init=cc.
- 2026-11-03 ~14:00: the tests were starved by the per-user CPU caps (norm 600, csbdevel 642) held by the RNase P and other jobs. With the user's approval, IgM workers move to `unlimited` (shares norm's nodes, own 432 cap, qos=normal). From 14:35, `igm_cryolo/rebalance.sh` runs every minute and moves pending IgM workers to whichever of norm, csbdevel and unlimited has spare cap and idle CPUs. It skips any partition where another job of the user is pending, so the RNA work keeps priority.

## 3D test outcomes and production (2026-11-03 evening)
- T2 (user's 15 classes, 192,481 particles, solve3D c1) won: final 11.5 Å (FSC 0.143), `~/Downloads/igm_T2_core15_final_pproc.mrc`. T1 (150k random of sel190) 15.6-16 Å. T3 (nstates=3) master OOM-killed at stage 5 (43.6 GB vs 16 GB request), dropped; stage-5 states 62,573/34,268/53,159 at ~16-17 Å in Downloads. T4 (EMD-13921 ref) cancelled at 21.8 Å by user.
- Core refinement 54937190: refine3D_auto on T2/1_solve3D/core15.simple, project volume (T2 recvol, 600 px), automsk=yes, regpass greedy nspace 5000 at 11.47 Å; run dir T2/3_refine3D_auto. Master 64 GB.
- Production `T2/prod.sbatch` (not submitted): after the add-on (2_solve3D_addon) ends, keep core + 2D classes whose add-on median corr >= core p25 (>= 30 scored particles), selection mkdir=yes, count check, refine3D_auto automsk=yes.
- First add-on (core + 200k mixed) OOM-killed at 43.6 GB before its report; its per-particle corr could not rank classes (junk scored above kept). Reran as two 64 GB add-ons, frozen = 183,255 updated core particles (9,226 never-updated core rows excluded):
  - drop (100k of the 110 dropped classes): IMPROVED, union 10.72 vs base 11.47 Å; cohort-only 12.34 Å, corr 0.74 with base. `~/Downloads/igm_addon_drop100k_only_pproc.mrc`.
  - keep (100k of the kept non-core classes): IMPROVED, union 10.72 Å; cohort-only 10.90 Å, corr 0.975. `~/Downloads/igm_addon_keep100k_only_pproc.mrc`.
- Core refine3D_auto: budget 32 iterations (~4 updates/particle, 25k sampled per iteration); iteration 1 9.48 Å, iterations 2-3 10.9 Å; pose-overlap convergence is bypassed under prob_neigh posterior draws. Master peak RSS 71 GB vs 64 GB request.

## Overnight plan 2026-11-03/04 (user: "work autonomously overnight... best IgM possible that I can plug into the tree"; code changes allowed)
- Core refinement 54937190 cancelled by user decision (stalled at 10.9 Å over iterations 2-8; refine3D_auto has no annealing: prob_beta_end is workflow-set and the parser rejects unlisted keys). Iteration-1 map 9.5 Å `~/Downloads/igm_core_refine_iter001_pproc.mrc`, iteration 8 `..._iter008_pproc.mrc`.
- Facts: T2 map centred (COM ~9 Å off), density radius p99 ~214 Å; EMD-13921 (600 px, 1.09 Å/px) density radius p99 ~287 Å, so it includes the Fabs. No raw movies (aligned DW micrographs only), so no RELION polishing. RELION 5.0.1 module exists; GPU queue mostly free. RELION track deferred: the tree needs a SIMPLE project.
- A1 54941349 `T2/neigh283k/`: refine3D refine=neigh nspace=20000 box_crop=432 filt_mode=nonuniform automsk=yes update_frac=1 maxits=10 on a copy of the keep add-on's published union (283,255 = core 183,255 + keep 100k). Test of deterministic local refinement.
- A2 54941350 `T2/addon_more/`: solve3D_addon chained on `addon_keep/t2keep.simple` (frozen 283,255) + 400k more kept non-core particles (seeded shuf) = 683,255, nparts=20, 96 GB.
- Next: final local refinement of the 683k union with the recipe A1 validates; then tree-readiness check (FSC stable, map, blind re-alignment control).
- A1 neigh283k: iteration FSC 0.143 8.84 -> 7.27 -> 6.23 Å (iter 3), overlap .34/.16/.71, mean pose change 19.5/23.8/16.5°. The FSC curve has a step at ~0.063 1/Å (phase-randomization onset) that puts FSC 0.5 at ~16 Å in iteration 3; the 0.10-0.20 1/Å band rose, halves independent (combine_eo=no). Master peak RSS 50 GB. Maps `~/Downloads/igm_neigh283k_iter001_pproc.mrc`, `..._iter003_pproc.mrc` (432 px, 1.514 Å/px). User: "looks good so far".
- A2 addon_more COMPLETED 1 h 51 min: IMPROVED, union 683,255 at 10.22 vs 10.72 Å (stage-8 cap), corr 0.995; published `T2/addon_more/t2more.simple`. Master 47 GB.
- A3 54942755 `T2/neigh683k/`: A1 recipe on the 683k copy, maxits=12, nparts=20.

## 2026-11-04 04:12 — A3 progress
- A3 (683,255, 54942755) passed its project check (all active rows scored). FSC 0.143: it1 9.34, it2 7.52, it3 **5.84 A** (overlap 0.24, 0.12, 0.68). Master peak RSS 49.9 GB.
- A1 (283k): it4 5.79 A, overlap 0.83, mean pose change 10.7 deg; still moving.
- Delivered ~/Downloads/igm_neigh683k_iter003_pproc.mrc (432 px, 1.514 A/px).
- Watcher bmuz8ktgn: fires at A3 iteration 7 or when A1 ends.

## 2026-11-04 04:25 — worker OOM; A3 resumed as A3b
- A1 (54941349) FAILED in iteration 5: worker part 7 OOM-killed (exit 137, 89 GB peak vs its 48 GB request, 3 workers on cn140). Last good A1 map: iter004, 5.79 A.
- Cause: at the 3.99 A matching limit, workers reach ~100 GiB, but job_memory_per_task=48000 let Slurm pack 3 per 256 GB node.
- A3's iteration 4 was at 3.99 A with 3 workers on each of four csbdevel nodes (shared with RNA). I cancelled it by ID: master 54942755 and workers 54944522-541.
- A3b (54944571, T2/neigh683k_r2.sbatch) resumes from A3's iteration-3 project (1_refine3D/t2more.simple, 04:09, copied to t2more_it3.simple). It sets job_memory_per_task=200000 (one worker per node) and runs maxits=9; output goes to T2/neigh683k/2_refine3D.

## 2026-11-04 05:17 — A3b iteration 1 clean
- A3b it1 (= A3 iteration 4): FSC 0.143 5.89 A (flat against 5.84), overlap 0.72, mean pose change 14.9 deg. Workers one per node at 200000M; peak RSS ~109 GB. No map delivered (no gain).
- It2 runs at the 3.99 A matching limit. Watcher bahornain fires after it3 or at job end.

## 2026-11-04 06:47 — A3b at 4.64 A
- A3b FSC 0.143: it1 5.89, it2 5.27, it3 **4.64 A**; overlap 0.72, 0.75, 0.86; mean pose change 14.9, 11.1, 7.7 deg. Matching limit 3.99 A (finest bank member).
- Delivered ~/Downloads/igm_neigh683k_r2_iter003_4.6A_pproc.mrc (432 px, 1.514 A/px).
- Watcher fires after it6 or at job end (maxits 9).

## 2026-11-04 09:02 — A3b at 4.22 A; final check staged
- A3b it4 4.48, it5 4.42, it6 **4.22 A** (overlap 0.876, mean pose change 4.1 deg). Delivered ~/Downloads/igm_neigh683k_r2_iter006_4.2A_pproc.mrc.
- The automask covers only the core (0.6% of the box, radius <= 140 A), so alignment is core-driven and the Fab arms average out. box_crop is a Fourier crop: same 654 A field of view.
- Staged T2/final_check.sbatch; watcher bxggx0bmf submits it when A3b COMPLETES:
  - (1) reconstruct3D of the final poses at 600 px, no automask or NU filter (Fab arms);
  - (2) blind control: refine=greedy on a 3% sample against the final reference, with angular change per particle from compare_poses.py.

## 10:45 — high-res branch staged (user approved)
- A3b iter 7: 4.09 A (overlap 0.885, 3.7 deg). Iter-7 map in Downloads.
- Matching handoff = finest NU bank member (static ladder ends at 4 A) or the ML pair at FSC if finer; so 4.0 A was the matching limit while FSC > 4 A.
- Change (separate build ~/SIMPLE_igm_hires, shared build untouched): lowpass_limits gains a 3.5 A rung; init_nu_filter drops trailing rungs that clamp to Nyquist (without it the evidence stage throws "not strictly coarse-to-fine" when Nyquist >= 3.5 A — caught by unit_reconstruction). Patch T2/nu_ladder_3p5.patch. Mac build: unit_reconstruction 19/19 PASS. Cluster build uses -march=skylake-avx512 (norm mixes x6150/x6342/x6740P).
- hires.sbatch: copy of A3b final project, refine=neigh box_crop=504 (1.297 A/px), maxits=5, job_memory_per_task=200000. Chain (Mac bg job be34n8o7s) submits it after A3b COMPLETED + final_check ended + HIRES BUILD OK.
- VPN/DNS to batch.ncifcrf.gov down at 10:20; chain waits for it.

## 11:05 — A3b COMPLETED (6 h 37 min)
- iter 8: 3.99 A; iter 9: 3.96 A, overlap 0.838 (down from 0.885), mean pose change 3.9 deg; not converged by SIMPLE's criterion.
- Final map: ~/Downloads/igm_neigh683k_final_iter009_3.96A_pproc.mrc (432 px, 1.514 A/px). Project for the tree: T2/neigh683k/2_refine3D/t2more_it3.simple.
- VPN back ~11:00 (Mac `host` lookup still fails while ssh works; chain no longer checks DNS). hires build job 54947297 submitted 11:04; final_check pending submission by watcher bxggx0bmf (10-min poll).

## 13:10 — final check COMPLETED; high-res submitted
- 600 px reconstruction of A3b's poses: FSC 0.143 = 3.96 A (same as the 432 crop). pproc in Downloads (igm_final_600px_3.96A_pproc.mrc).
- Blind control (refine=greedy, update_frac=0.03, 20,498 particles re-searched globally vs the final map): projection direction median 0.00 deg, 96.8% within 3 deg, 97.6% within 6, 1.5% beyond 30. Full rotation 96.4% / 97.3% / 1.7%.
- Tail (>30 deg) histogram + shared-axis test (pseudo-5-fold flips) submitted: job 54948024 (tail_hist.py, tested on synthetic 72-deg flips).
- High-res refinement submitted by the chain: job 54948021 (box_crop 504, 3.5 A rung, separate build).

## 13:25 — VERDICT: poses tree-ready
- Tail >30 deg: 340 of 20,498 (1.66%), angles spread 30-180 following the random-rotation density (rising to 180), axis concentration 0.40 (isotropic 0.33), no 72/144 peaks -> no pseudo-5-fold flips; tail = junk/low-signal, left for the tree.
- Tree input: T2/neigh683k/2_refine3D/t2more_it3.simple (683,255 ptcls, 3.96 A).

## 13:30-19:40 — IgM hetero tree launched, root OOM, resumed
- Tree dir /mnt/beegfs/elmlund/afan/igm_cryolo/tree/igm683k_s1 (launcher igm_tree.sbatch: production f692_lf keys, mskdiam=600, norm x6342, nodepar 8, job_memory_per_task=64000 on the project copy).
- build_igm = math_reworks/src_lf + tree_stage_box.patch: stage crop boxes grow with field of view so Nyquist <= lp/1.25 (IgM 110/164 px; RNase P keeps 64/88 exactly). Fast gate 15/15.
- Root: mean, flex (8 states, 33-143k; maps in Downloads igm_tree_root_flex_state1-8, states differ at 150-300 A radius), cheap, fflex (mask 92 % of support), stA0, recB0, site discovery (2 sites) OK in ~5.5 h.
- 19:16 site batch OOM: one-job steps sized 8 GB + 0.025 MB/ptcl (RNase P box 160) = 25 GB; IgM used 67 GB.
- build_igm2 = + tree_mem_rate.patch (SIMPLE_TREE_MB_PER_PTCL env, default 0.025 unchanged). Resume job 54949975 (igm_tree_resume.sbatch, MB_PER_PTCL=0.12 -> 90 GB at the root): root set back to step 18 (site batch), OOM'd dirs moved to *.oom1; site jobs resubmitted at 89991 MB.
- hires (54948021): iter 1 3.94 A at matching 3.50 A; iterations ~3 h under the shared CPU cap.

## 2026-11-05 02:21 — IgM tree root split
Root split: r1 = stage-B states 1+3 (588,416, IgM; core-shift difference merged), r2 = state 2 (94,839, junk blob). Depth 1 started on r1 (8 flex states, focus mask 45% of support). Stage-B maps in Downloads (igm_tree_root_stageB_state01-03). hires iter3 3.85 A delivered; iter 4 running.

## 2026-11-05 03:23 — hires iteration 5 FSC inflated; continuation cancelled
Iter 5 reported 3.50 A but its FSC sits ~0.2 out to 3.52 A and falls off a cliff exactly at the 3.497 A matching limit (halves share information to the matching limit); 4-5 A FSC dropped (0.675 -> 0.655), cFAR 0.44 -> 0.39, overlap 0.93 -> 0.87. Iter 4 (3.74 A, smooth tail) is the best SIMPLE result. hires2 (54952235) cancelled after 1 min, no output dir. RELION CtfRefine chain (relion/rchain_igm.sbatch) goes on the iter-4 poses (relion/probe dump); pilot gate 54952185 running. Iter-4 maps were overwritten by SIMPLE; iter 3 (3.85 A) and iter 5 maps exist.

## 2026-11-05 04:37 — RELION baseline
Full 683k RELION reconstruction from SIMPLE iter-4 poses: frame cc 0.999, white; masked phase-randomised FSC 3.85 A, smooth tail (no cliff) = honest number for the best SIMPLE poses (SIMPLE own FSC said 3.74). ctf_refine running in relion/full (job 54952302).

## RELION CtfRefine (job 54952302, done 06:54)
- Full CtfRefine (per-particle defocus, per-micrograph astig, beam tilt, odd n<=3, even n<=4 per 39 day-groups): 3.94 A vs 3.85 A without — it HURT.
- Defocus moves are sane (sd 172 A, 4/683k > 2000 A); B ~0, scale ~1. Even-aberration fits scatter widely day to day (4th-order terms -34..+26) -> suspected noise.
- Ablation running: 54953401 (A: refined defocus/astig only), 54953402 (B: A + beam tilt + odd). Best raw map goes to Downloads.
- Downloads/igm_relion_ctfref_3.94A_raw.mrc is the WORSE map; the 3.85 A uncorrected sum is full/Cons/sum_w.mrc.
- Ablation: A (defocus/astig only) and B (+tilt, odd) both 3.940 A, same as full CtfRefine -> the refined defocus/astig itself costs ~0.04 FSC from 5 A out; tilt/aberrations neutral. SIMPLE's per-particle CTF beats RELION's refit.
- BEST RELION map: Downloads/igm_relion_cons_3.85A_raw.mrc (= full/Cons/sum_w.mrc, iter-4 SIMPLE poses, original CTF).

## Tree depth 1 (08:22)
- r2 (94,839, junk side of the root) -> leaf r2_g1_2: the tree did not split the junk further.
- r1 (588,416 IgM) stage B done: 3 states, 70 % of the cross-half variance in 10.6 % of the molecule (localized = real-split signature). Maps: Downloads/igm_tree_r1_stageB_state01-03.mrc. Split test 54953899 running.
- r1 split test -> LEAF r1_g1_2_3 (588,416). EM groups 328k / 134k / 127k; difference power over half-map noise P/N 1.68 (1|2), 1.26 (1|3), 2.22 (2|3); 1|3 and then 1|2 under the 2.0 effect floor -> merged. No quality split (|d| < 1), no pose-flip flag.
- Tree result: two leaves, IgM 588k (no discrete state above the effect floor) + junk 95k. IgM variability is real but small/continuous at this depth.

## 2026-11-05 ~10:30 — why the IgM tree misses the Fab arms; workflow folder
- Workflow folder: $B/WORKFLOW (README + copies of scripts/refs/class lists + link_* to projects/maps).
- r1 candidate sites all core (<= 84 A). r1 stage-B maps are ZERO beyond ~100 A: the node focus mask (falloff_mask) keeps vs > MFRAC(0.2) * max(vs); the dense core sets max.
- r1 mask_variance per 20 A shell: p99/noise 340-720 in the core, 40-70 at 100-140 A, 9-19 at 140-200 A, 5-7 at 200-280 A, all < 0.08 of max -> mask 0 beyond 100 A. Integrated over shell volume, more reproducible variance lies beyond 100 A than inside. Root mask was whole sphere (junk vs IgM variance everywhere).
- Cheap test: build_igm3 (= build_igm2 + SIMPLE_TREE_MASK_FRAC override in falloff_mask only, default 0.2; tree/tree_mask_frac.patch), job 54954335. Then igm_mfrac.sbatch x2: 100k seeded subset of the r1 leaf, maxdepth=0 site_layer=no: TAG=stock vs TAG=f002 (MASKFRAC=0.02). Dirs tree/mfrac_stock, tree/mfrac_f002.
- 100k box-600 tests (54954361/2) cancelled by ID (+ workers 54954363-78): box 600 too slow (user: test isolated steps; box size is the issue).
- Isolated mask step offline (mask_frac_shells.py on r1 mask_variance, dilation 3 vox): share of shell kept, f=0.2/0.05/0.02/0.01: 100-120 A 0/8/31/60 %; 140-160 A 0/0/3/18 %; 180-200 A 0/0/0/3 %; >200 A ~0 at all f. A lower MFRAC can NOT reach the arms (arm vs ~0.005 of max, near noise): only a mask covering the arm envelope by construction can.
- Test bed tb100k (54954389): same 100k r1 draw, selection + prune_project + scale_project newbox=200 (3.27 A/px). Next: root split at box 200, stock vs arm-envelope mask.
- tb100k DONE in ~12 min (prune + scale_project newbox=200; 100,000 active, poses kept, 33,878 stacks in tree/stack_parts_sc): tree/tb100k/igm_sub_sc.simple.
- tb_tree.sbatch (maxdepth=0, site_layer=no) on the test bed: tb_stock 54954620 (build_igm3, stock mask).
- build_igm4 (54954631) = build_igm3 + SIMPLE_TREE_MASK_FILE (supplied focus mask replaces falloff_mask's; tree/tree_mask_file.patch). Masks from the RELION consensus at box 200 (mk_tbmasks.py, 54954632): tbmask_env (envelope incl. arms), tbmask_arms (env outside r 100 A). Watcher submits tb_env and tb_arms once gated.
- 683k IgM tree 54949975 COMPLETED 11:01 (56,082 s): 2 units (IgM 588,416 + junk 94,839). No site: r1_s1 pooled to r, 0 poles pass halves there; r1_s4 refined ranking accepts no pole; r1_s2/s3 no pole. Final maps in igm683k_s1/1_hetero_tree/final_maps.
- tb masked runs: tb_env 54954647, tb_arms 54954648 (build_igm4); masks: envelope to ~220 A, 37-60 % of the 100-200 A shells.

## 2026-11-05 ~12:30 — pentamer vs hexamer (user: "I only get hexamers")
- Ring-centred angular power of maps projected down the disc normal (symfold.py + local centring): EMD-13921 n5 0.60/0.76, n6 0.00/0.04 at 90/110 A (clean 5-fold + an elongated rod = subunit-1 F(ab')2 held by the J chain). Ours: solve3D core15 n6 0.71/0.79; RELION 3.85 n5 .34/.23 n6 .40/.61 -> 6-fold, no rod. DJ 15_solve3D (12.1 A) is a rod, not a disc. 2D stars: 92 5-fold, 138 mixed, 246 6-leaning.
- Paper (Chen et al. 2022 Nat Commun, PMC9588798): myeloma IgM pentamer, SEC removed hexamers (<5 % of IgM); crYOLO 1.82M -> 2D 546,583 -> 5 ab initio + heterogeneous refinement, 2 best classes 301,286 -> NU-refine C1 4.3 A; 3DVA 5 conformations + subtraction/focused refinement on subunit 1. Asymmetry from the J chain.
- Reading: our map = pentamers averaged over ~60-degree registers of the pseudo-hexagonal core (core-only automask cannot see the J chain / Fab-1), not real hexamers. The blind control tested 72/144 flips only, not 60. Tree results on these poses are compromised.
- tb3d_reg (refine=greedy update_frac 0.2 vs EMD-13921 at box 200, lp 10, automsk=no): INCONCLUSIVE. Frame-change cluster 0.4 %; residual rotations follow the random-rotation density (rising to 180), axes isotropic (0.34). greedy refine3D does not re-search globally from a foreign frame (and the hand may differ), so no register can be read. A real test needs dock_volpair (both hands) first. Decision deferred to the de novo full-mask runs tb3d_s1 / tb3d_s3; symfold_tb3d.sbatch queued by a watcher on their completion. Note: tb_3d.sbatch reg branch reads 2_refine3D (1_ is the scale step).
- Box-200 root-split tests DONE (1 h 38 - 1 h 53 each vs 15 h at box 600): all three -> leaf. stock: mask 100 % of support at this root, 70 % of stage-B variance in 35.5 %, no site split; env: groups 57k/43k, P/N 1.115 < 2 -> merged; arms: 53k/47k, P/N 1.087 -> merged (variance in 3.3 %). No quality split (|d| < 0.11). Masking the arms in does NOT produce a discrete split at 100k on the current (possibly misregistered) poses.
- 15:xx EARLY READ of the full-mask de novo runs (box-200 test bed, automsk=no, random 100k of the 683k): tb3d_s1 stage 3 (~27 A) n5 0.60 vs n6 0.12 at 100 A = PENTAMER. tb3d_s3 stage 2: states 1/2 (11k/12k sampled) n5 0.31/0.42, n6 0.02; state 3 (6k) 6-fold at 120 A (n6 0.39). So the 6-fold consensus came from the core15 de novo start (side-view-heavy 15 classes), not from the data. Figure scratchpad igm_pent_vs_hex_early.png. Final maps + fold check on watcher byfgt97qu.
- PLAN (user agreed): pentamer and hexamer analysed separately; particle count data-driven. pool200 job 54957909: copy of igm_all.simple (all 2,400,639, no 2D preselection) -> dupcheck.py (pairs per micrograph < 50/100/150/200 A) -> scale_project newbox=200 into igm_cryolo/pool200/stacks_sc. Next: pentamer/hexamer/junk 3D classification on pool200 seeded by the full-mask de novo maps; then ResLog (random subsets reconstructed at box 600 from fixed poses) per population to set its size.
- FINAL de novo full-mask (ring-centred analysis; the COM-centred mid-run read was misleading): 1 state (14.9 A) drifted back to 6-fold (n6 0.73 at 110 A). 3 states: s1 43.7k no ring symmetry (side/junk), s2 29.9k 5-fold (n5 .36-.39) WITH an elongated rod like EMD-13921's Fab-1, s3 26.4k 5-fold (n5 .53 at 110 A). No 6-fold class among the 3. Reading: a single model averages pentamers over pseudo-hexagonal registers -> 6-fold; classification recovers the asymmetric pentamer. Hexamers not yet shown to exist in this data. Figure igm_denovo_final_rings.png. pool200 rerun 54960476 (nparts=60 after worker time-limit failure, fail1 kept).

## 2026-11-05 ~16:10 — sentinel mode (user: monitor all IgM chats, work autonomously for hours)
- pool200 DONE (rerun 54960476, nparts=60): pool200/igm_pool_sc.simple, 2,400,639 ptcls, box 200, 3.27 A/px, 358 GB stacks; all 2.4M active in ptcl2D and ptcl3D.
- 5-state de novo on all 2.4M (user-approved): job 54961035, pool200_s5/ (solve3D nstates=5 automsk=no mskdiam=560 c1, nparts 20 x nthr 16, 48 GB/part). Watcher runs symfold_s5.sbatch on completion.
- Sentinel: scratchpad sentinel.sh (new FAILED/TIMEOUT/OOM/NODE_FAIL under igm_cryolo, 10-min checks, 4 h life).
- Other threads: DJ session (dj_ws2): 2D of 2.25M gave no star class; re-classifying 612k disc particles (54961009) + pick-overlap check (54961013). Fixture session: fixture build 54960473.
- Memory igm-11077-consensus-poses-tree-ready marked SUPERSEDED (6-fold register average).

## 2026-11-05 ~18:10 — sentinel check
- s5 (54961035, pool200_s5): iteration 23 of 133 after 65 min (solve3D samples 10k per state = 50k per iteration, 8 stages MAXITS 20/20/17/15/12/12/12/25, so every particle is visited ~2.8 times). Sampled populations at iteration 22: 10.9k/6.6k/10.8k/12.6k/9.1k; only state 4 has a real FSC so far (cFAR .39, the rest .01-.05). Estimated end 6-9 h. Watcher bmt8v958t runs symfold_s5 at the end.
- DJ session: pick overlap done (DJ under-picks the star views ~2x; discs are 17 % of OUR 2D too, 138k of the old 683k consensus came from disc classes). Disc re-classification 54961009 running. Its 10 x 32-cpu workers were on csbdevel (disc2d.sbatch's update_project lacked qsys_partition=norm); nothing was queued on csbdevel (1,980 idle cpus, no RNase P jobs), so it was left to finish and the session was told to set norm from now on.
- Fixture session: fixture built (realism matched within .02); 4 current-tree runs + real_root_stock on norm. Its note's "6-fold class ~20-26 %" came from my wrong mid-run read; corrected by message (final ring-centred 3-state read: no 6-fold class).
- No new IgM job failures since the sentinel started.
- 18:15 sentinel fired: fixture (c) subtraction jobs 54963622/23 FAILED in 10 s (exit 7; no cons_chunk_*.mrcs, a script bug). The fixture session had already relaunched them as 54963668/69, now past that step (rest map built). Sentinel re-armed (4 h).

## 2026-11-05 19:30 — user: investigate autonomously, picking adequacy not established, csbdevel cores available
- 5-state s5 correction: independent multi-state solve3D = 5 stages (MAXITS 20/20/17/15/12 = 84 max), lpstop 6.5; stage 2 converged at iteration 25. Since stage 3 (box crop up) ~10 min/iteration: ~2-4 min in workers (one part lags ~4 min), ~3.5 min serial master assembly (40 s per state). Expected end ~00:00-03:00.
- PICKING AUDIT (audit/, csbdevel) on DJ's 1,000-micrograph test batch (= micrographs 1-1000 of igm_tpick.simple, >= 10 sessions):
  - pickers: ours (full-data template picks, pickrefs_v3), seg (picker=segdiam moldiam_max=420, ~300 blobs/mic; union keeps diameter >= 150 A), gau (picker=new moldiam=300), t3d (192 refs: 32 directions x 3 in-plane of the de novo pentamer+rod map tb3d_s3 state 2 and the 6-fold map tb3d_s1; refs3d.mrcs, 400 A mask, 20 A), dj (DJ's workspace_2 extraction).
  - jobs: refs 54965606, seg 54965519, gau 54965520, t3d 54965690; watcher submits pa2.sbatch (union within 50 A -> extract box 600 on a copy of DJ's ds_batch project -> box 200 -> solve2D 150 classes -> pa_score.py class IgM score (best correlation vs 1,200 projections of 3 maps, in-plane search), view tilt, ring n5/n6 -> per-picker recall/junk + overlays).
  - pa_full2d (54965691): same class scoring on our full-data 300-class 2D (all 2.4M picks) = IgM share of all picks, and what the disc classes are.
- REPLICATES (reps/, csbdevel): rep_A 54965521, rep_B 54965522 = disjoint seeded 300k draws of the 2.4M box-200 pool, solve3D nstates=6 full mask lpstop=10, symfold at the end.
- DJ session told to run its 3D the same way (random 300k of DJ's 2.25M, 6 states, lpstop 10, csbdevel).
- ~20:00 FOLD READOUT CORRECTION: ring n5/n6 at 110 A flips with the centre (EMD-13921: centre .13/.15, COM .05/.14, min-n1 .25/.42, max-n5n6 .74/.01). Centring-free |F| angular test (audit/fold_ft.py): EMD-13921 45 A n10 .32 / n6 .01 (pentamer; core n12 .53 at 35 A); our consensus maps + tb3d_s1 n6/n12 dominant, n10 ~0 (6-fold average stands); tb3d_s3 states 2/3 no clean n10 -> their "5-fold" call is unconfirmed.
- SIMPLE symmetry_test (audit/symcal, lp 20, msk 400, C2-C8): not discriminative (EMD-13921 c5 .83 = c6 .83 -> c1; all maps C4-C7 within .02, C2 top). Not a yardstick.
- full-data class scoring v1 non-discriminative (all classes .67-.98: radial profile dominates); v2 removes ring means (angular structure only), calibrated on core15/keep175/drop110/disc: job 54966126.
- ~20:30 full-data class scoring v2 (angular structure only): core15 median .92, keep175 .65, drop110 .42, disc .37. Classes >= .6 = 52.8 % of all 2.4M picks; of all picks: side views 43.3 %, tilted 8.4 %, face-on 1.1 %. Discs score as non-IgM (but may be face-on IgM blurred by in-plane ambiguity; untested).
- Star classes 92/138/246 (the only face-on views): cleanly 6-FOLD by every measure (6 outer blobs; |F| n12 .88-.95 at 80 A, n10 ~0; both centrings agree). Consistent with hexamers OR with pentamers whose in-plane rotation locked onto the pseudo-hexagonal core (5 Fc + J chain) and averaged over 60 deg, in 2D as in 3D. 6-fold appearance cannot decide.
- PAPER (Chen 2022, same EMPIAR-11077): commercial myeloma IgM pentamer, SEC removed "trace amounts of dimers and hexamers"; crYOLO 1,819,931 picks -> 2D 546,583 -> 5 ab initio + heterogeneous refinement, 2 best classes 301,286 -> NU-refine 4.32 A (symmetry not stated); J chain locally 3.5 A; 3DVA on subunit-1 F(ab')2 hinge; no hexamer class, no orientation/2D-class description, no discussion of the pseudo-hexagonal core.
- DECISIVE TEST = J-chain slot: jreg (job 54967433, jreg/): 150k random particles from IgM-scoring classes, box 200, refine3D_auto vol1=EMD-13921 (b200) ref_pose_init=cc, full 560 A mask, automsk=no, C1. If the J chain + Fab-1 rod come back in both halves beyond 15 A -> registers resolved -> then classify the sixth slot with fixed poses (J chain vs a sixth Fc).
- 21:05 5-state 2.4M run CANCELLED (user): master 54961035 + 15 workers by ID; watcher stopped.
- Disc test (audit/full2d_radial.npz, audit/disc_profiles.png): rotational-average profiles of the 51 disc classes lie between side views and face-on stars: compact (to ~90-100 A radius), smooth, NO 120-180 A Cmu2/Fab halo that the face-on stars (92/138/246) and the EMD-13921 face-on projection show. Consistent with face-on Fc cores with disordered outer domains OR junk of that size; 2D cannot decide. Plan: cross-tab disc particles vs replicate 3D states/poses.
- EMD-13921 at box 200 down its disc normal: pentagonal core + the subunit-1 Fab rod (audit/emd_proj_pair.png). Sixth-slot plan changed to an unsupervised 2-3-state classification from the jreg poses (no synthetic hexamer reference).
- 21:36 REGISTER TEST jreg DONE (54967433; 150k IgM-class particles, box 200, refine3D_auto vol1=EMD-13921 ref_pose_init=cc, automsk=no, 25 its, not converged, overlap .37): FSC0.143 15.6 A, cFAR .41-.44. Rotational CC with EMD-13921 about its disc normal (20 A): final/even/odd .886/.884/.885 at 0 deg vs .63 at 60/72/120 deg; EMD self 1.0 vs .65; C6-averaged EMD control .83 flat at 0/60/120/180. => both half-maps carry EMD-13921's asymmetric register (J chain + Fab-1 rod side), not a 6-fold average. Caveat: started from EMD-13921 and stayed near 15 A, so reference bias not excluded.
- Control jreg_c6 (54972007, jreg/3_refine3D_auto): same particles/settings, vol1 = EMD-13921 averaged over six 60-deg rotations (emd_c6_b200.mrc, no register info). If the asymmetry returns, the particles carry the register.
- 22:15 USER approved the per-particle 60-deg register fix of the 683k (after the C6-start control confirms). Implementation staged in regfix/ (regfix.sbatch, REF=<registered map>): dock the registered pentamer (both hands, regdock.py picks by 20 A corr) onto the 683k consensus (T2/neigh683k/2_refine3D/recvol_state01.mrc = box 432, 654 A -> 1.514 A/px), scale to 600, then refine3D refine=neigh (its coarse pass searches all coarse directions + full in-plane, so a 60-deg register can switch) nspace=20000 box_crop=432 filt_mode=nonuniform update_frac=1 automsk=no lp=12 maxits=3 from the consensus poses. Register flips measured with tree/register.py (before/after ptcl3D).
- Particle count (user: "are there more or less than 683k"): IgM-scoring 2D classes hold 1.27M (>= .6) / 0.80M (>= .8) picks; the 683k includes 138k disc-class particles of unknown identity. Plan: after the 683k fix, align the remaining IgM-scoring particles + discs against the registered pentamer and keep those that fit -> data-driven count, then ResLog.
- DJ session told the register fix is implemented here (no duplicate).
- 22:25 DJ session: dock3_ref.mrc + dock3.json in ph_tests/ (reference's pseudo-6-fold axis on the consensus axis; MIRROR hand; slot masses Fc 70-91, J slot 9). CORRECTION: its face-on "no hexamers" verdict is weaker: the synthetic hexamer was built about a 5.4 px-off axis; a correctly built hexamer vs C6(pentamer) = .965 (not .90), so the real .983 only leans to pentamers. Rerunning tonight. Register-fix margins give a further hexamer handle: a hexamer (no J chain) scores all six registers equally -> near-zero margin in clear face-on/tilted views.
- Control jreg_c6 at iteration 12 of 25 (~23:20 end); interim symmetry-breaking check on its iteration-12 map and halves running locally.
- 22:35 control jreg_c6 at iteration 12: STILL SYMMETRIC (best CC per 60-deg sector .698-.723, top-2nd .009, map and both halves) vs jreg .878 vs <= .61 (top-2nd .269). Design flaw: a perfectly C6 reference scores all six registers equally per particle, so nothing drives symmetry breaking; the control cannot confirm the user's condition. Reference-free evidence of register information = DJ's fine face-on 2D (9.4k gapped pentamer particles + Fab-1 rod).
- 100k register-fix TEST submitted (54974230, regfix/tb): REF=ph_tests/dock3_ref.mrc AXIS=dock3.json; checks: j uniformity, margins, gold-standard halves, J chain/rod beyond the reference detail. The 683k fix waits for the user on this result.
- 22:50 DJ session RETRACTED the face-on "no hexamers" verdict: corrected synthetic face-on 2D (54974178, hexamer H2 on the true axis, P:H 2:1) purity .744 vs chance .662; six-lobed classes are P/H mixtures, a hexamer-rich class looks identical. Our 7.9k six-lobed face-on classes are UNDETERMINED. Next (approved here): a 7th candidate (hexamer H2, full subunit wedge incl. Fabs) in the register scorer, validated on a synthetic 3D mix (P:H 2:1 and 1:0 control) with confusion by view and margin; if it validates, add it to the real tb100k/683k scoring.
- regfix_tb first attempt 54974230 FAILED (old register_fix.py without CHUNK args); resubmitted with regfix_handover/register_fix.py as 54974704 (tb.fail1 moved aside).
- 22:39 REGISTER FIX TEST ON REAL tb100k FAILED (54974704, regfix/tb): j 0..5 = 20664/28854/20025/4626/6548/19283 (not uniform); FSC0.143 orig 7.43 A -> corr 12.34 A, corrtop 12.82 A; |F| corr MORE 6-fold (80 A n12 .44 vs .19; 45 A n10 .01 vs .15). Diagnosis job rotuni (54976452): all particles turned by the same 60 / 180 deg with the scorer's rewrite; unchanged 7.4 A = axis/rewrite fine, per-particle choice is the problem; a drop = axis/frame wrong. DJ session told. The 683k fix is NOT run.
- Picking audit step 2: union 161,573 (ours 64.6k, t3d 61.4k, seg>=150A 61.3k, gau 35.3k, dj 41.5k on 650 mics); ours found 76 % of t3d's picks, 51 % of DJ's; 130,852 extracted; scale failed (scale_project workers race to mkdir dir_target) -> resumed as pa2b (54976025) with stacks_sc pre-created; IgM cutoffs recalibrated to .5/.6/.7.
- 23:00 RETRACTION: the "6-fold register average" verdict on the consensus is WRONG. Self-corr of the tb100k consensus (rec_orig) rotated 60/120/180 deg about the dock3 axis: 30-8 A .34/.26/.21, 12-7 A .04/.01/.01; after the wrong-hand fix .70/.62/.61 and .21/.17/.17. The consensus is mostly ONE register; the 6/12-fold |F| = pseudo-hexagonal core. The DJ session also found EMD-13921 AS-IS (not mirror) fits rec_orig better (.483 vs .421 at 30-8 A), so the dock3 (mirror) reference made the register scores meaningless. The rotuni control (7.43 / 7.52 A) validates only the pose-rewrite math, not the axis (any rigid rotation keeps the FSC) - my design error.
- Retry (54977294, regfix/tb3), input-free: score each particle against the six rotations of the consensus ITSELF, cross-validated (eo=0 vs odd half, eo=1 vs even half); corr + corrconf (change only where margin >= median); FSC + self-corr + |F|.
- DJ session 7-way synthetic (matched models, real views, 3 deg pose noise): P->P .973 (false-H .027), H recall .820, register accuracy .953, overall .922 (chance .664); needs a matched reference on real data.
- Memory corrected (igm-11077-consensus...).
- 23:05 DJ synthetic de novo (54968346): pure-pentamer 1-state stays pentameric (ring .68/.02 max-centred; |F| n10 45 A .12 vs model .29 = partial register mixing, like our rec_orig .15). P:H 2:1 2-state solve3D: NO separation (state 1 = 44,997, purity = chance) even with the easier v1 hexamer -> a missing 6-fold class in rep_A/B (or DJ's runs) is only WEAK evidence against hexamers. The separating test = per-particle 7-way scoring with a matched reference. Asked DJ whether P2_cons/H2_cons are built from our consensus (candidate 7 for real data).
- Control jreg_c6 final: 25.2 A, cFAR .156 (never broke symmetry; moot given the retraction).
- 22:59 REGISTER VERDICT (regfix/tb3, 54977294): cross-validated scoring of tb100k against the six rotations of the consensus itself: j=0..5 = 97549/937/317/160/204/833 (margin quartiles .98/1.42/1.81); confident half: 49964 j=0 vs 38 other. FSC orig 7.43 A, all calls 7.60 A, confident-only 7.43 A; self-corr unchanged (.34 at 30-8 A). => THE CONSENSUS IS REGISTERED (>= ~97.5 %); no register fix for the 683k. The "6-fold" signals were the pseudo-hexagonal core + the centring-sensitive ring readout. Pentamer vs hexamer stays open -> 7-way with a matched hexamer candidate (asked DJ about H2_cons).
- 23:15 CAVEAT on the tb3 verdict (DJ): rec_orig wedge masses about the axis 43.9/32.5/29.0/23.9/35.8/63.9k = smooth one-sided pattern (EMD: five equal Fc + empty J slot) -> either real (J chain + Fab-1 side) or the axis POINT is 1-2 px off; an off point penalises every j != 0 and could manufacture the 97.5 % j=0. The hcons refinement changed only the direction (1.15 deg), not the point. Checks: (1) axis scan (54977840, regfix/axscan): point +-2 px along two perpendicular directions + refined direction; (2) asked DJ for an axis-free face-on check (average the six-lobed vs gapped face-on sets with CONSENSUS poses).
- Matched hexamer candidate from OUR consensus: ph_tests/hcons/H_orig_cons.mrc (C6 of one rec_orig subunit wedge to r 230 A incl. Fabs; box 200, consensus frame). DJ's synthetic P2_cons/H2_cons are EMD-derived (mirror, no Fabs) = not matched.
- 23:20 DJ: equalising the five Fc slots would need a ~20 A centre move (contradicted by the centred hub) -> the uneven slot masses are real one-sided density; the axis-independent support for "registered" = max 60-deg self-corr over the axis search .44 (rec_orig) vs .82 (scrambled). Matched hexamer on exactly the dock3 axis: ph_tests/hcons/H_orig_dock3axis.mrc. Delegated to DJ (it owns the 7-way scorer + calibration): real-data 7-way on tb100k, cross-validated halves, D distribution vs synthetic, H calls at |D| > .33/.63, and the decisive check = reconstruct the called-H subset (binary patch, not import_particles) for six equal Fc slots and no J gap in both halves.
- 23:30 AXIS SCAN (54977840): keep-register dock3 .9755 (38 confident switches), refined dir .9751 (39), point +-2 px along u/v .9936-.9954 (0-3 confident switches). Moving the point off makes switches RARER (off-axis penalises j != 0); dock3 is the best point tested; no hidden minority anywhere. With the axis-free max 60-deg self-corr .44 vs .82 -> the 683k consensus is registered (>= ~97.5 %); no register fix.
- 23:45 PENTAMER vs HEXAMER ANSWER for the consensus set (DJ, job 54978042, ph_tests/real7/): cross-validated 7-way on tb100k: D quantiles 5/25/50/75/95 = -3.31/-2.15/-1.49/-.96/-.35 (no H mode); called H 1.2 % (D>0), .22 % (D>.33), .04 % (D>.63) = at/below synthetic pure-P false-H; by view face-on 4.0 % (vs 2.5 % false), tilted .8 %, side 1.1 %. Most-H-like 2,000 vs strongest-P 2,000 reconstructions: same lopsided pentamer (slot masses/mean subH 1.13/.80/.84/.73/.86/1.64 vs subP 1.11/.87/.76/.72/.91/1.62; 60-deg self-corr .565 vs .569; halves agree). => hexamers < ~1 % in the consensus. Registers re-confirmed (97,791 j=0). Caveats: consensus particles only; face-on excess too small to reconstruct; face-on pose-averaging (54978096) running.
- 00:00 USER: the 3.85 A map "looks like a hexamer". CORRECTION of the 23:30-23:45 verdicts: J slot / Fc slot in our consensus .59-.70 vs ~.11 in EMD-13921 -> J chain smeared over the ring = core register scrambled (or hexamers). tb3 ("97.5 % keep register") and real7 ("no hexamers") used references derived from the consensus itself = self-referential; they cannot detect a register average baked into the reference. Solid: 2D face-on gapped pentamers (9.4k); with consensus poses the gap vanishes (DJ). Decisive route: J-gap-bearing reference, correct hand, J/gap/rod mask, success = J slot empties in independent half maps; DJ's reference-free per-class face-on angle test (54978217) running.
- 00:25 USER OK to start the full alignment without asking. J-slot metric (audit/jslot.py, EMD frame, ring r 30-110 A, |h| < 40 A, six 60-deg sectors phased on EMD's J sector): EMD-13921 J/Fc .24 (.22-.23 at 30/20/15/10 A low-pass: the gap is a LARGE feature, so a 30 A reference still imposes it); first jreg run (full EMD ref) .42 in full and both halves; C6 controls 1.00-1.11.
- So the reference must NOT carry the gap: emd_jfill_b200.mrc = EMD-13921 as-is with its J sector (r < 118 A) filled by the neighbouring Fc rotated 60 deg (J/Fc 1.15), Fab-1 rod (r > 118 A) kept; low-passed 30 A. A J gap re-emerging in independent half maps = data-driven register.
- Cancelled (by ID) the gapped-reference runs: jreg_lp30 54978326 (+ workers) and fullreg 54978327 (+ 40 prune workers); dirs moved to fullreg.cancel1 and jreg/4_refine3D_auto.lp30_cancelled. Resubmitted with the J-filled reference: jreg_jfill 54978458 (150k validation, ~1 h) and fullreg 54978459 (all 1,643,183 candidates = 2D classes >= .6 + 51 disc classes; prune + refine3D_auto CC init, nsample 150k; ~8-12 h).
- 00:35 DJ axis-free face-on check (reference-free): consensus-pose averages of the gapped face-on set lose the gap (n10 .04 vs 2D .25-.40); consensus in-plane angle minus 2D gap angle per gapped class: ~68 % main register, ~32 % off (mostly +-60). 7-way calls by view: j != 0 face-on 22.8 % (low margins), tilted 1.5 %, side 1.3 %. DJ reads side/tilted as registered. NOT established: those side/tilted calls are self-referential, and if ~92 % were J-registered the consensus J slot would be ~.12 of an Fc slot (DJ scale), not .59-.70. Decisive = J-filled-reference runs (54978458 / 54978459).
- 00:40 PICKING AUDIT RESULT (pa2b 54976025; audit/audit_score.tsv, audit_overlay.png): 130,852 union particles extracted on 1,000 micrographs; IgM classes (score >= .6) 54 classes, 48,040 particles. Raw recall of all IgM: ours 51.2 %, t3d 48.9 %, dj 33.9 % (650 mics), seg 31.5 %, gau 18.5 %; junk share ours 55 %. "Missed by ours" 23,425 = mostly off-centre duplicates (our nearest pick < 150 A for 59 %; > 250 A only 5,492 = 11.4 % of all IgM) -> corrected recall ~89 % (~85 % at 200 A). Face-on IgM 3.7 % of all IgM for every picker (sample orientation, not picker bias); ours clearly misses ~27 % of face-on vs ~11 % overall. Verdict: picking adequate, no full repick; optional union with t3d picks (+~10 %).
- 00:50 DJ reference-free face-on registers (54978407, ph_tests/facecons/facereg_registers.txt, corrected t2more_it3 poses): 32 % of 9,235 gapped face-on change register; pooled J/Fc .22 -> .11, |F| n10 .04 -> .14; gap angle 198.4 = consensus lightest slot 197.7 = EMD as-is J 195.8 -> the consensus has the J slot in the right place, only partly emptied.
- DJ 3D fix with EMD as-is clean-gap reference (54978511) FAILED: calls dominated by 120-180 deg; 12.3-15.2 A; J/Fc 1.10 (filled), self-corr .85. Cause: EMD's heavy slot is the Fab-1 rod root at J-60; ours is heavy at J+120 and light at J-60 -> the score matched mass, not the gap. My jfill reference keeps EMD's rod, so it may inherit the pull; my readout watches all six sectors (any sector emptying in both halves = data-driven register). If none empties: cancel fullreg, rerun with DJ's gap-only reference (C6 ring of our consensus minus the J wedge, r < 125 A) after its validation on the 9,235 reference-free face-on registers.
- 01:00 DJ gap-only reference (ph_tests/gaponly/G_gaponly.mrc, J/Fc .18): face-on validation agreement .783 (chance .17; all-j0 .678; margin >= q50 .908). On tb100k as a scorer: 68 % change, FSC 9.9-13.1 A, J/Fc .43-.55, self-corr up to .78-.85 -> side/tilted register calls unreliable with every scorer tried; only face-on gaps are readable.
- jfill 150k (54978458) at iteration 3 (23 A): J sector .92 (ref 1.12), halves equal; J+180 sector .82 (unchanged). Pending end (~1 h).
- Proposed to DJ: hexamer test on the 7,958 six-lobed face-on particles: G-scoring margins vs the gapped set (misaligned pentamers find a gap; hexamers give near-zero margins), then average with confident calls (gap -> pentamers; six equal lobes -> hexamers).
- 01:15 USER: "is this a shortcoming in the 2D? ... don't want code changes, a correct workflow". Answer: mainly 3D (core-only automask hid the rod/Fab cues, leaving the small J chain to set the register); 2D only shows it for face-on views and fine 2D already resolves most. Workflow test (no code change, input-free): selfreg (54979140, regfix/selfreg): tb100k from consensus poses + consensus map, refine3D refine=neigh (coarse pass = all in-plane angles) nspace=20000 filt_mode=nonuniform update_frac=1 FULL 560 A mask automsk=no, 8 iterations; judged by the J sector emptying in both halves (regfix/jslot_cons.py, consensus frame on dock3 axis: rec_orig .56 (DJ .59), scrambled 1.19).
- 01:20 SIX-LOBED HEXAMER TEST (DJ, 54979324, ph_tests/sixlobe/): the 7,667 six-lobed face-on particles are MOSTLY GAPPED PENTAMERS: gap excess .66-.82 vs the gapped set (C6 no-gap controls ~0); cross-fitted averages (register from one Fourier half, placement in the other) J/Fc six-lobed .24-.40 vs gapped .09-.21 (no gap = 1). Upper bound: <= ~25 % of six-lobed (~1,900 = ~11 % of 16.9k face-on) show no gap = hexamers OR gap-invisible pentamers. The six-lobed 2D classes are not hexamer evidence. DJ corrected its 2D window numbers (gapped face-on .16 -> .08 before/after register fix; EMD .13, consensus .60 on its 3D windows).
- 01:20 csbdevel at the account cap (642): CANCELLED fullreg 54978459 (+18 workers, by ID) to unblock the decisive tests (jfill 150k, selfreg); its pruned 1.64M project in fullreg/ is kept for a restart with the validated recipe. USER: "you can shift some jobs elsewhere, norm or unlimited" -> moved 12 pending workers (selfreg + DJ pool200 r300k) to `unlimited` (scontrol update); all running.
- 01:55 WORKFLOW TEST FAILED (selfreg 54979281): tb100k from consensus poses + map, refine=neigh, FULL mask, automsk=no, lp 10, 8 its: J/Fc .56 -> .60 (both halves), FSC 7.43 -> 10.2 A, cFAR .12, overlap .99 (poses stuck). Full-mask refinement does not self-register and costs resolution (flexible Fabs). jfill 150k (54978458) at iteration 21: 24 A, cFAR .23, J/Fc .87 (ref 1.12), no other sector empties -> the cancelled fullreg recipe would not have worked.
- Summary of register methods tonight: self-referential scoring (circular), EMD mirror / as-is scoring (rod mismatch -> scramble), gap-only scorer (reads face-on, scrambles side), full-mask refinement (stuck), J-filled refinement (24 A, no gap), reference-free face-on (works, face-on ~4 % of IgM). Side/tilted views do not give per-particle register with any method tried.
- Next route (paper-like, no code): class-based selection - multi-state de novo (rep_A/rep_B 6 states, full mask) -> keep classes whose map shows a clear J gap. Frame-free readout regfix/jslot_self.py (2D core-slab projection, centre = best 60-deg ring repeat): EMD .13, C6 1.00, consensus .63, jreg .41.
- 02:10 REPLICATES DONE (rep_A/B: 300k random of the 2.4M each, solve3D 6 states, full mask, lpstop 10; ~16 A). Frame-free J gap (regfix/jslot_self.py: lightest sector / Fc; EMD .13, consensus .63, C6 1.0):
  - rep_A: s1 15.2k .48 (ring self-corr .32, edge centre) | s2 59.1k .91 (.85, symmetric) | s3 25.3k .68 | s4 10.2k .58 | s5 145.4k .28 (.78; |F| n10 .16-.20 at 45/35 A) | s6 44.8k .55
  - rep_B: s1 44.6k .89 (.89, symmetric) | s2 117.3k .26 (.60) | s3 12.1k .31 (ring .24) | s4 31.9k .75 | s5 16.5k .17 (ring .38) | s6 77.6k .59
  => each replicate has ONE large class with a clear J gap (A5 48 %, B2 39 %) next to symmetric register-mixed classes (A2, B1 ~.9): the paper's class-based route separates a registered pentamer.
- Class refinements (unlimited partition): classref_A5 54980810, classref_B2 54980811 = refine3D_auto state=S from class poses + class map, full mask, automsk=no; judged by the gap in half maps at higher res. Then dock A5 vs B2 for reproducibility.
- 02:57 Class refinements: classref_A5 (54980810) 18.7 A, J readout .25 in map and both halves; classref_B2 (54980811) 15.6 A, .19 (ring centre at the search edge). High-res local refinement of A5 submitted (classhr_A5 54982348: prob_neigh_mode=geom pose_cont=yes automsk=yes).
- 03:09 2D candidates vs de novo classes (reps/clsoverlap): fullreg/states.txt (1,643,183) keeps 89.6 % of A5 and 93.9 % of B2, 8-21 % of the junk-like classes (A3, B3-B5). Four disjoint ~308k candidate draws were staged (reps/canddraw.sbatch) but NOT submitted (see 03:20).
- 03:15 READOUT TRAP: jslot_self.py took the disc normal from second moments of the whole thresholded map, so the arms set it (maps docked into one frame came out 45-50 deg apart). Core-only normal (r < 125 A, regfix/jdir2.py): EMD .12 (clean single gap), EMD-C6 1.00, J-filled .83, consensus .63, A5 .38 (gradient), B2 .19 (gradient, centre at the search edge); small junk classes A1/B3 also read low.
- 03:20 RETRACTION of the 02:10 class-based "breakthrough". Deposited-frame test (reps/regemd, regfix/regemd.py): every map docked into EMD-13921 through a 280 A core mask (20 A) and read on EMD's own axis/centre; register contrast = core CC(k=0) - best CC(k!=0). Calibration: EMD .189, EMD J-filled .075, EMD-C6 .001. Ours: consensus -.020, A5 +.013 (mirror +.005), B2 -.061 (mirror +.032), all 12 replicate classes -.034..+.038, i.e. at the C6 level. In the full-mask docking into A5, EMD-C6 fits better than EMD (.864 vs .818). Caveat: our cores match EMD worse than EMD-C6 does (core CC .73-.79 vs .895), so register averaging vs a structural difference near the core is not separated. classhr_A5 cancelled by ID (master + 10 workers; it had slipped 13.6 -> 21 A).
- Where this leaves IgM: no 3D map of ours (consensus or any de novo class) shows the deposited J gap. DJ's synthetic tests say the register IS recoverable per particle with a matched reference (.95, side views included), so the missing piece is a matched registered reference. Candidate routes (not started): (a) the paper's own recipe in RELION (Chen 2022: 5 ab initio classes + heterogeneous refinement, keep 2) on rep_A's 300k; (b) C6 symmetry expansion + one-sector subtraction + focused classification without alignment on a 100k consensus subset, validated first on DJ's synthetic set.

## 07:35 (2026-11-06) — s6 half-map check verified; independent replicates launched; DJ tests s6 as the matched scorer
- Verified DJ's s6 half-map output (54983521): even +.138, odd +.137, full +.138; k=0 best in both halves. solve3D aligned both halves to one shared merged NU reference (l_lpset), so this rules out reconstruction noise, not reference imprint.
- Independent check (running): two disjoint random 300k of DJ's box-200 pool, rows 300,001-600k (`reps/djrep/rB`, 54984558, csbdevel) and 600,001-900k (`rC`, 54984559, norm) of the same seeded shuf. Identical solve3D (6 states, 560 A mask, lpstop 10, nparts 10 x nthr 24), then regemd on every state in both hands. ~9 h each.
  - A J-gapped state (contrast > +.075, k=0 best) in either draw = s6 confirmed from disjoint particles and independent de novo references, and the route scales: each draw is also a production chunk of registered particles.
  - Our own pool gave 0 of 2 draws with a register state (rep_A/B); DJ's gave 1 of 1.
- DJ session asked to test s6 as the matched reference in its per-particle scorer on tb100k.
  - Deciding readouts: agreement with the 9,235 reference-free face-on registers (>= .85, vs gap-only .783), and FSC not worse than ~7.6 A (orig 7.43 A).
  - Supporting readout: J/Fc in both halves.
  - No 683k fix or refine3D until read.
- Sentinel: no IgM job failures in the 4 h to 07:20.

## 08:40 (2026-11-06) — RETRACTION: s6 is not the registered pentamer; regemd is a > 40 A readout
- 08:20 sentinel: DJ's s6 scorer job 54985371 (ph_tests/s6score2) stopped at its own docking gates.
  - EMD control vs P_asis .716 (gate .90).
  - s6's own axis 5.9 deg / 9.2 px off EMD's.
  - s6 vs EMD .06 at 30-8 A; s6 vs rec_orig .149, against EMD's .464.
  - Docked directly into the consensus, s6's light sector sits 158 deg from the consensus J; EMD's sits +1 deg.
- Band check (job 54986056, regfix/regband.py; core r < 130 A, CC per band vs EMD at k*60 and vs EMD-C6):
  - 60-30 A, k=0: s6 .26 < C6-level s1 .30 (consensus .51, EMD-C6 .77).
  - s6's relative 60-deg drop is .62; EMD's own is .37, and a noisy or register-mixed EMD copy cannot exceed .37. So s6's asymmetry is lopsidedness beyond the J register, which docking lays on EMD's J side.
  - 30-20 A: all 18 de novo class maps (12 ours, 6 DJ's) are ~0 vs consensus .40. solve3D's "16 A" is shared-reference imprint; DJ's shell readout agrees (real resolution ~40-60 A).
- TRAP: regemd's 20 A-Gaussian contrast reads > 40 A structure. A register claim needs three things:
  - (i) a relative 60-deg drop <= EMD's;
  - (ii) EMD-like CC at 30-20 A;
  - (iii) docking into the consensus putting the J at the consensus J window.
- s6 retired as the register reference; DJ session told to stand down on s6.
- The reps/djrep replicates (54984558 csbdevel, 54984559 norm; rB already had worker 54986133) lost their premise. My scancel by ID was blocked by the auto-mode classifier, so they run until the user decides.
- Route left (user's call): RELION C6 symmetry expansion + focused classification of one slot on the consensus poses (no code change), validated first on tb100k against the 9,235 reference-free face-on registers. Fallback: the paper's de novo recipe; our de novo maps argue against it.
- 08:50 CORRECTION (DJ session): "J 158 deg off, fit .149 vs EMD .464" came from DJ's v2 placement (54985371), which failed its own gates; it is not a measurement. Criterion (iii) is untested for s6, not failed. The verdict rests on the band check alone: s6's 60-30 A relative drop is .62 vs EMD's .37, and it has no 30-20 A structure. DJ's v3 shell CCs agree: vs EMD .43/.21/.09/.06 at 60-40/40-30/30-25/25-20 A; halves .82-.90 down to 9 A = imprint.

## 09:40 (2026-11-06) — USER approved the RELION symmetry-expansion pilot; job 54987705 (igm_cryolo/symexp/)
- Set: tb100k (regfix/tb/rows.tsv, consensus poses at box 200) + the 9,235 gapped face-on particles with reference-free registers (DJ facereg; consensus poses from cons_ptcl3D /3, face2d stacks), about 109k. Duplicates between the two are allowed.
- Poses moved to the dock3 ring-axis frame: R' = R Q^T, s' = s - (R c)_xy (SIMPLE convention; SIMPLE e1/e2/e3/x/y = RELION rot/tilt/psi/origin, per export_relion).
- Six copies per particle, R' Rz_k. Copy k equals DJ register j = k (checked: Q^T Rz_k Q = axis_rot(n, 60k)).
- Steps, in one job:
  - self-test;
  - build;
  - relion_preprocess to box 100 (6.54 A);
  - expand;
  - relion_reconstruct of copy 0 and copy 1;
  - gate G1: RELION map vs rec_orig in the frame >= .8 at 30-15 A, ring axis on z, copy-1 map = Rz_1-turned map, J slot agrees;
  - slot mask: soft 60-deg wedge at the J slot, r 25-115 A, |z| < 45 A;
  - relion_refine_mpi K=2 --skip_align, tau2_fudge 4, 25 iterations (16 ranks x 6 threads, csbdevel);
  - score: register call = the copy with the highest J-class probability, agreement with the face-on registers, calls by view, particles with no J copy;
  - box-200 half maps, consensus vs called registers: FSC0.143, J/Fc in both halves.
- Success: face-on agreement >= .85 (gap-only scorer .783, all-j0 baseline ~.68); FSC not worse; the J slot empties in both halves.
- 10:12 pilot G1 PASS (job 54987757; earlier attempts 54987705 and 54987717 stopped on fixable steps):
  - RELION 4 relion_preprocess no longer operates on STAR files, so a Python Fourier crop to box 100 replaced it.
  - The consensus "lightest slot" is not a usable J slot (masses .69/.65/.69 nearly equal), so the J slot now comes from the face-on registered map.
  - RELION map from the re-framed poses vs rec_orig: .996 (60-30 A) / .973 (30-15 A). Ring axis on z (60-deg self-corr .584 vs .033 about the box axis). Copy-1 map = Rz_1-turned map (.993, next .546).
  - Face-on particles with reference-free registers: one empty slot at 185.5 deg, J/Fc .20, the five others 1.07-1.24 of the mean. With consensus poses: J/Fc .41.
  - Consensus slots from that J onward: .65 .69 1.11 1.93 .92 .69 (heavy slot opposite the J).
  - Mask: 60-deg wedge at 185.5 deg. Classification running from 10:12.

## 11:00 (2026-11-06) — RELION symmetry-expansion pilot FAILED on all three criteria (54987757 + 54987915; symexp/)
- K=2 skip-align classification in the J-slot wedge (25 iterations, ~25 min; no signal subtraction) split the copies by the consensus HEAVY LOBE, not J vs Fc.
  - Minority class: .191 of copies. It matches the consensus turned to copy 2 (CC .89 in the mask). Both class maps have the same slot density (1.27 vs 1.33 sd units); neither looks like an empty J (CC .21 / .26 vs the face-on J copy).
  - Mask placement was right: in the face-on registered map, the J copy has the lowest density in the mask (-.12 vs .25-.40).
  - The heavy lobe sits ~147 deg from the J (frame ~332-337 deg). It is tied to the consensus frame, not to each particle's J.
- Face-on agreement:
  - J-class rule: .008, with call minus truth at +2 (.35) or +3 (.28).
  - Lobe-angle rule (offset fitted): .386; resultant length .555; calls 0 or 5 whatever the truth.
  - All-j0 baseline .678; gap-only scorer .783.
- Half maps (box 200, RELION): consensus poses FSC0.143 7.11 A, J/Fc .61/.60 in the halves, 60-deg self-corr .344. Called registers: 7.88 A, J/Fc .64/.63, self-corr .574. Worse, the same signature as every earlier scorer.
- Assets kept for a revised run: pilot_z.star (box 200), pilot_z_b100.mrcs/.star, expanded_b100.star, faceon_j/0_b100 maps, slotmask_b100.mrc, cls2/.
- Candidate next step (user's call): relion_particle_subtract (keep the slot) + K=3 classification seeded J / Fc / heavy from face-on subset 1, validated on face-on subset 2.

## 12:10 (2026-11-06) — revised RELION pilot FAILED for side/tilted views (54988649, symexp/)
- Pipeline:
  - one-iteration consensus model;
  - relion_particle_subtract of all 655,410 copies, keeping the J-slot wedge, recentred on (-76, -7, 0) A, box 40 at 6.54 A;
  - K=3 skip-align classification seeded J (4,623 face-on subset-1 copies at k=j), Fc (9,246 at k=j+-1) and heavy lobe (20k tb at k=2); tau2_fudge 4, 25 iterations, ~10 min.
- No J class formed:
  - class 1 (J seed): .270 of copies, CC with seeds J/F/H .56/.71/.66;
  - class 2: .508, .55/.75/.74;
  - class 3: .221, .99 with the heavy-lobe seed;
  - centre density 1.141 / 1.131 / 1.262 sd, none empty.
- Face-on agreement:
  - held-out subset 2: .785 (all-j0 .678, gap-only .783); seeding subset 1: .788;
  - by margin quartile: .58 / .73 / .87 / .96;
  - by truth: j=1 .73, j=5 .69, j=2-4 .19-.27.
- tb tilted/side: median margin 0.00; >1 J copy for .91-.98 of particles = ties, no calls.
- Half maps with the calls: FSC0.143 7.11 -> 8.84 A, 60-deg self-corr .344 -> .737 (scrambled), J/Fc .61 -> .45 per half.
- VERDICT: fixed-pose classification (with or without subtraction) reads the register in face-on views only (~.79, like the gap-only scorer) and not in side/tilted views. Together with every scorer tried and the de novo runs, per-particle register is not recoverable this way here.
- The paper's route differs in kind: multi-class refinement WITH full re-alignment (cryoSPARC ab initio + heterogeneous refinement, 4.3 A).

## 12:20 (2026-11-06) — USER DECISION: the IgM tree input is the current consensus AS IS, with the register caveat
- Input: tpick_all/3D_tests/T2/neigh683k/2_refine3D/t2more_it3.simple (683,255 particles; FSC 3.96 A; RELION 3.85 A). The earlier tree run tree/igm683k_s1 used exactly this input and is the IgM baseline.
- Caveat to state with every IgM tree result:
  - about 1/3 of particles are in the wrong 60-deg register (face-on 2D: 32 % off, mostly +-60);
  - the consensus J slot is ~60 % filled (J/Fc .61);
  - the consensus is locked onto a heavy lobe ~150 deg from the J;
  - a split by J register or a "J-slot state" is this nuisance, not biology;
  - hexamers < ~1 %.
- Register work closed: fixed-pose register reading fails for side/tilted views (RELION pilots v1/v2, every scorer). The paper-style multi-class refinement with re-alignment remains an optional longer run.
- The reps/djrep replicates (54984558/9) serve no purpose; cancelling them needs the user's OK.
- 12:35 USER: cancel the djrep repeats. Cancelled by job ID: masters 54984558 (rB, csbdevel) and 54984559 (rC, norm), then rB's workers 54989115-54989124. The queue is empty after; partial outputs stay in reps/djrep/rB and rC (not read).
- 14:45 sentinel: the 4 h watch expired with no IgM failures, but sacct is refusing connections (slurmdbd on fsitgl-head01p:6819), so its check was blind for the outage; logs checked instead. Fixture site-split gates: gate 0 (54996933) NORMAL STOP 14:32; the gate-1 no-crop rerun (54999269) finished 14:33 in 5 min: r1 .86 (was .78), z1 13.9, agreement with Python .82 (was .80, needs >= .9), still no call vs Python not discrete. The fixture session got its wake-up at 14:34 but has not acted on it. Nothing IgM is queued, so the sentinel is not re-armed.
