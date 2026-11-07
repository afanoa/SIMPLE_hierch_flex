---
name: pcg-focus-mode-and-vol-keys-2026-10-20
description: "pcg_focus=yes (fixed-background focused PCG, x = x_fix + P delta, x_fix = per-half consensus of all states) and vol1..volMAXS + MAX_CMDARGS 200 implemented 2026-10-20 in worktree vigilant-cartwright-dcd646 (uncommitted); batch2 build ~/SIMPLE_focus, snapshot snap/focus1; test arms on fam6"
metadata:
  node_type: memory
  type: project
  originSessionId: 3ca8858f-79b1-4b46-8b2f-59410f52ac90
  modified: 2026-10-21T03:42:22.056Z
---

Uncommitted diff in the Mac worktree `.claude/worktrees/vigilant-cartwright-dcd646` (base origin/master 5703662e7).
- `pcg_focus=yes` (refine3D, reconstruct3D; needs rec_backend=pcg + pcg_mskfile as the focus mask; refuses trail_rec,
  automsk, envfsc, pcg_solvent, atm_prior, nonuniform, state=). Master (simple_rec3D_pcg_strategy) first solves the
  consensus pair from ALL states' raw accumulators per half (add_raw_accum_weighted, sphere support; ml pair under ml_reg;
  written as recvol_focus_bg.mrc), then each state solves P H P delta = P(b - H_data x_fix) (reconstructor_pcg
  set_fixed_background; window_mul/div add/subtract x_fix; data-only kernel apply). State FSC + replay P_tau = FSC of the
  deviation x - x_fix. A state lacking a half ships the consensus.
- vol keys: simple_args_generator.pl emits vol1..volMAXS (99, parsed from simple_defs.f90); simple_cmdline MAX_CMDARGS
  100 -> 200 (K=24 worker line had 74 keys, so >50 states would have hit 100) + argcnt guards.
- batch2: worktree ~/SIMPLE_focus (detached 5703662e7 + copied files), build ~/SIMPLE_focus/build (BUILD_TESTS=ON,
  build_focus.sh N -> makeN.log/buildN.done), snapshot ~/SIMPLE_focus/snap/focus1 used by toy55/fam_classify_focus.sbatch
  and fast_comp_focus.sbatch (copies; FOCUS=yes, PCG=yes, NTHR). Swept swing mask: rnasep_klem/ff3_swing30hs/focus/swing_swept.mrc
  (toy55/swing_mask.py; contains 98 % of region_5, 6x larger). Invariant check: toy55/focus_check.sbatch (RD, MASK).
- RESULTS (fam6 of ff3_swing30hs, 8 bins seeds, refine3D prob_neigh geom lp 10, nspace 5000 box_crop 88, snap/focus2, toy55
  fam_classify_focus3.sbatch EXTRA): focus + swept mask 5.2..25.8 A (sd 3.2-4.1, min state 966) | gridding control 6.6..22.3 |
  unmasked PCG 7.5..23.8 (ends drained to 39/390) | focus + region_5 12.1..16.7 (partial mask) | hard region_5 13.6..15.0 (one
  state empty). refine3D_states baseline (5e1e6aeee) 6.9..22.7. Invariant: native box 160 states == consensus bit for bit outside
  the mask; at box_crop 88 the Fourier-cropped mask's ringing tail leaks ~1e-5 rms (read_and_crop, pre-existing, also widens hard support).
- K=32: toy55/fast_comp_focus.sbatch KSEED=32 (snap/focus1, norm -c 48) -> ff2_swing30hs/comp_k32 NORMAL STOP, 6 its, 2 h 3 min
  (44 min of it = initial NU filtering of 32 states, 65 s each), purity 0.60 seeds -> 0.97, smallest state 30. K=24 had died on vol21.
- MEMORY TRAP: refine3D (plain) defaults to NSPACE 20000; prob_tab_neigh holds pfts_refs to Nyquist (interpklim) + ft_ref/ft_ref2 for
  nstates*nspace refs: 8 states x 20000 at box 160 = 76 GB PER PART (4 local parts = 300 GB); fam_classify's 120G only "works" on
  csbdevel 2 TB nodes (not enforced), norm 256 GB nodes OOM. nspace 5000 + box_crop 88 -> 7.6 GB per part.
- One focus bug found on data: l_solvent_weight read before initialization in the consensus ml pass (fixed, build 4).
- Build 5 (snap/focus3, toy55/fam_classify_focus4.sbatch): pcg_focus floors the read_and_crop'd focus mask at -min(mask) (the
  resampling ringing depth; 0.021 at box 88, 0 uncropped). Swept rerun: 4.2..23.9 A, sd 3.3-3.7, min state 2186; states now differ
  in 2.0 % of the box (was 48 %), 99 % of it inside the mask, nothing more than 5 voxels from it.
- TRAP for PCG unit fixtures: analytic projections FFT'd with SIMPLE's normalized FFT reconstruct at BOX x the density
  (2D FFT /N^2 vs 3D /N^3); generate data with pcgop%forward_plane of x*get_invenv() to stay on the operator's scale.
Related: [[klem-bottom-up-organization-2026-10-20]], [[never-push-to-origin-without-permission]].
