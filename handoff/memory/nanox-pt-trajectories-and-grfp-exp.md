---
name: nanox-pt-trajectories-and-grfp-exp
description: "Map of /mnt/beegfs/elmlund/NanoX (READ ONLY): the 15 published Pt LC-TEM trajectories, raw frames, per-frame poses; GRFP tests live in afan/grfp_tseries_exp (job 54727429, 2026-10-30)"
metadata:
  node_type: memory
  type: reference
  originSessionId: 36dc059e-25bd-4e47-9f7c-87d5f98d1a8b
  modified: 2026-11-01T06:11:41.123Z
---

User 2026-10-30: NanoX is READ ONLY ("read only. nothing else"); experiments go in a new dir in
/mnt/beegfs/elmlund/afan (here: grfp_tseries_exp, tiers 1-3 of the GRFP test plan, sbatch 54727429 on norm).

**The 15 trajectories** (NanoX/time_series.txt; = the 15 of Meana-Paneda 2025 Sci Adv): 71-17 ptcl1-6,
71-6 ptcl1-4 (raw full frames exist: full_frame_data_sets/71-17/mrc, 71-6/frame*), 4D NP87 and
NatComChem np3-6 (tracked stacks only). 0.358 A/px, box 160, 400 fps (frame names count per second).
- per particle: full_frame_data_sets/<s>/NP<s>/3_tseries_track_particles/tracked_ptcl_N/ holds raw,
  subtr, den stacks, .box (one row per frame) and reproc/<n>_autorefine3D_nano/ptcls_oris.txt
  (indstk = frame, per-frame e1-e3, x, y, corr); take the highest-numbered refine dir.
- segment tables (frame ranges + #ATOMS) are in reproc/log, else parse_seg_inf.csv + segvols/*_ATMS.pdb.
- TRAP: reproc imported DENOISED stacks (subtr_den_cong / *_den, kPCA across frames): their poses and CCs
  carry denoiser-induced inter-frame correlation; temporal residual tests must use the raw stacks.
- TRAP: filetab.txt paths point to old /home/elmlundho/data/... locations; remap by basename.
- raw frames: 1920x1792 float32 with integer values (~70-480/px/frame) -> gain unknown until tier 3.
- TRAP: in some reproc runs (71-17 ptcl4-6) indstk=0 and the frame number is in pind.
- Results 2026-10-30 (t1_prop1, t2_textstats_v2, t3_rawframes, t1_mra in the exp dir):
  poses (denoised, all 15): P(within 5 deg) at lag 1 = 0.3-0.8 (71-17, NatCom), 0.01-0.09 (71-6, NP87);
  e-fold ~10-20 frames; jumps near-uniform; 230-3300 distinct views/movie -> no long stay-or-jump runs.
  S = CC^2/(1-CC^2) 0.05-0.53 (median 0.11); a 0.5 s event timing needs ~3-38 atoms (median 11);
  atom loss 0.3-9 /s. Raw frames: integer ADU, var/mean ~28 per electron, ~6 e/px/frame (~1.8e4 e/A^2/s),
  neighbour-pixel corr 0.30-0.34 (per-pixel binomial thinning invalid), background lag-1 acf 0.07-0.11
  (white = -0.005), still 0.01-0.04 at lag 10. Prop 1 smoother 1-1/sqrt(1+2u) confirmed (0.782 at u=10).
  MRA toy: iid sigma exponent 6.2; oracle ell exponent ~-1.6; fixed windows 4-14x worse than oracle
  when ell << sigma^2 (gap grows with sigma) -> windowed half of Conjecture 1 fails.
- Others' active work (ask first): NiPt/, CdSe/ (meanapanedar2). Simulated null: cdse_wurtzite_simulated
  (abTEM, rigid CdSe, Brownian rotation, 10k frames, 0.2645 A). Do not use heterogeneousNPs recs (info.txt).
- T1 final (t1_residuals_r2, job 54736659): raw-frame fit CC 0.10 (S ~0.01 full band); particle-only lagged
  excess 0-0.045 at lag 1, e-fold 9-35 frames (pose-linked), none in 71-17 p2/p4 -> no slow continuous
  fluctuation on Pt. Denoised stack inflates it 2.4x. Shift sign in these oris = "+x,+y" for every trajectory.
  TRAP: SIMPLE reproject needs an EVEN nspace even with oritab (build_refspiral) -> pass nspace=2.
- Dissolving_particles: NP01.mrc is Au51 (gold, 0.3169 A), not NP87; NP87 atom counts differ by pipeline (548 vs 346).
- User 2026-10-30: Au51 is NOT useful (never published) -> dropped; pose smoothing tests on 71-17 p3.
  Prefer published trajectories for pilot work.
- PILOT (user-approved scope 2026-10-30): 71-17 p1, p3, NP87; plan in
  doc/implementation_notes/grfp_tseries_pilot_plan_2026_10_30.md; no jobs until the user says go.
- Stage 0 done 2026-10-30 (pilot/stage0): NP87 raw frame i = NP87_subtr frame i+2000; fission gate FAILED
  (shared bg 5-12% after fitted fission; 50-frame blocks worse, 24-36%) -> improve the noise model, keep
  fission. Stage 0b (two time scales, block b-1 -> b): max|xcorr| < 0.01 all three; half-sum shared 4x4
  NP87 0.010 PASS, 71-17 p1/p3 0.05 (slow bg drifts between blocks); px level 0.012-0.015. In-band S (1-3 A,
  atom support) p1 0.039, p3 0.111, NP87 0.062.
- Stage 1 done (pilot/stage1): plain windows do NOT beat published segments (half agreement ~linear in frames;
  w passing 0.8 x segment agreement: p1 2000, p3 none, NP87 1000 -> ~6/3/8 maps vs 5/3/7 published). Old
  detect_atoms counts unusable on half maps; rerun with ~/SIMPLE_master1030 (origin/master 6ddc4a0d2) in
  pilot/stage1_atoms. TRAPs: import_particles needs mkdir=no to update projfile; swap_stack/reimport keep indstk.
  TRAP: S = CC^2/(1-CC^2) is a PER-PIXEL SNR, not the whole-particle SNR; never plug it into
  N/(S k) timing formulas. TRAP: the |P|-profile support finder clamps at 66 px on these reprojections.
- Stage 1b (new detect_atoms, pilot/stage1_atoms): early segments half-consistent (p1 376/375, p3 709/712);
  windows < 2000 frames and late segments noise-limited (counts inflate, |A-B| 50-160); only NP87 w2000 passes.
- Time-specificity (pilot/stage1_control): 1/2-2/3 of in-particle half agreement is time-specific, decays
  smoothly with gap -> gradual one-way change (dissolution). Outer shell: same-time-only bg, 2-3x larger p1/p3.
- Flex trajectory_denoise on halves (pilot/stage3_flex): CV picks h = 8 frames, tau 16, rank 32 (= max), neff ~11.
- NEW MODEL (user OK'd drafting 2026-10-30): doc/implementation_notes/grfp_tseries_dissolution_model_2026_10_30.md:
  per-atom absorbing occupancy (one detachment time T_a), mean-field EM, timing rms 5.1/s per atom (dm_cpconst),
  prototype scripts/dm_occupancy.py -> pilot/dm (probe pilot/dm_probe). NP87 Stage 1 halves = 8,338 posed frames
  (of 15,971 recorded, stride 1.92 -> convert to recorded frames before quoting seconds; dm_report.py does).
  Conventions verified: PDB/smpd = voxel index; project R x onto (x,y) = SIMPLE reproject (CC .977); frame pos =
  centred - shift. RESULT: toy passes all 3; NP87 real 128 atoms, median |T_A-T_B| 0.77 s (null 7.3 s), N(t)
  tracks 7 published segments, T vs coordination rho .80; p1 0.96 s but over-calls late; p3 only 45 atoms.
  Real s_eff 2-6x below toy from kappa 3-5 (pose-linked correlated misfit) -> next: pose refinement vs atoms.
  Flex halves: 2-3x half agreement at equal w (NP87 w250 .065 -> .213), but outer-shell shared bg up 2-4x too.
- 2026-10-31 (note 6f/6g): pose refinement HURTS real data (dropped). Static per-atom offsets reproduce A/B 0.9
  (~0.2 A rms) and give best timing (NP87 0.49 s); kappa is not positions. p1's ~150 "core losses" = signal fading,
  not departures: per-block contrast (scripts/dm_blur.py, pilot/dm_blur) removes them (148/160 -> 8/11), toy-checked.
  Contrast swings 2-3x block to block, shared by halves; costs NP87 timing (0.77 -> 2.1 s). Open: view-dependent
  contrast. Half agreement can be inflated by shared artefacts: also check order vs published segments.
- 2026-10-31 night (note 6h-6j): per-frame contrast real (halves 0.8 per frame; view explains ~25-40%).
  Fixed-contrast timing rides shared contrast drops (77% of NP87 departures) -> its 0.77 s was an artefact; toys
  with the measured contrast: modelling contrast cuts true timing error 2-5x. Departures are GRADUAL fades in
  place (scripts/dm_ramp.py, pilot/dm_ramp3: NP87 ~6,400 frames ~30 s; midpoint dT 0.60 s), confirmed by the
  published segment maps (dm_segcheck.py: ramp occupancy vs map density r .80/.69/.51, graded bins). Hazard follows
  CURRENT bonded coordination on NP87 beyond time and second-shell controls (dm_physics.py). p3 departures not
  confirmed by maps. Hazard prior (dm_ramp.py --hazard) recovers toy beta 0.7 but is NOT adopted: no toy or map
  gain; its 5x real half agreement is prior-driven. Model of record = ramp + per-frame contrast, uniform prior.
  TRAP: half agreement misled twice (fixed contrast, hazard prior) -> score timing on matched-SNR toys / maps.
- 2026-10-31 day (note 6m): NP87 = paper trajectory EFC2 (Meana-Paneda 2025 Sci Adv, 337 atoms; 4D_NP87 in
  time_series.txt). Wave model via SIMPLE master simulate_scattering (dm_wave.py, pilot/dm_wave): multislice ==
  projection for these particles; one global defocus explains less than the 2 deg view predictor. TRAP:
  simulate_scattering treats dfx=0 in an oritab as absent (draws 2 um) -> use |df| >= 1e-4 um. Pose errors exist
  (neighbour-pose weights 0.5 on NP87, far views 0) but correcting them leaves the contrast spread unchanged. The rest is a
  shared state persisting ~25-50 frames, tied to lateral motion; per-window defocus fits it (defocus drift or blur).
- 2026-10-31 evening (note 6n, dm_dfblur.py, pilot/dm_dfblur2): Python strong-phase imaging == simulate_scattering
  (cc 0.99999; Kirkland signs, Cs -0.011 mm). Physical image beats the Gaussian template held-out (+4-5/7-8/14-19%
  NP87/p1/p3); per-window df/B shape changes real (match drift toy on p1/p3); fitted blur tracks rotation speed on
  p3 (r 0.5-0.6). Contrast LEVEL still unexplained: df/B add nothing beyond a free per-window scale. TRAP:
  per-window df/B without a window scale acts as a contrast knob (misspec toy) - always include a misspec control.
- 2026-11-01 (note 6o, dm_ramp_phys.py, pilot/dm_phys2 + dm_phys2t): physical image IN the occupancy fit. Toys
  (physical truth) time fades 2.3-6.5x better than Gaussian (p3 Gaussian 39% false); real maps r -.015/+.009/-.036 and
  more unconfirmed core departures (NP87 34-37 vs 6) -> model of record unchanged. Defocus/sign family degenerate.
  TRAP: patch-s1 sweeps cost ~12-25 min at 32 workers (4x the probe); main sweeps ~30 s -> use >= 12 sweeps.
- NiPt survey 2026-11-01 (read only): STEM-HAADF graphene-liquid-cell PtNi from Jungwon Park's lab (SNU deck
  NanoX/NiPt/260911_PtNi_reconstruction_share_v3.pptx). Zips: Dissolution_#1-3 = DENOISED tifs (2,000/2,000/3,000
  frames; #3 has monomer_ground_truth_session.json), particle#1-3 = mrc (3,000/6,000/10,611). Ruben's static nano
  pipeline runs (analysis2D_nano -> autorefine3D_nano -> detect_atoms Pt,Ni) in /mnt/beegfs/elmlund/ruben/NiPt/<set>/.
  Pt/Ni species split not admissible except dis3; no time-resolved work yet. Others' active work: ask before using.
Related: [[grfp-2026-application]].
