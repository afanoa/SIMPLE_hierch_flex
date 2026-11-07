---
name: f692e-landscape-pruning-test-2026-10-20
description: "user 2026-10-20: f692e landscape (closed / partial / swinging arm, overlay image) is the target, but it oversplit and some volumes lacked signal -> masked gold-standard pruning of every f692e region board and zoom (splitbench8 --focus) + exact f692e reproduction f692r"
metadata:
  type: project
---

User verdicts 2026-10-20: the gold-standard trees (gs1/gs2/polish) lack diversity; the Oct 18 f692e landscape (overlay
of state maps: arm sampled over its swing, closed and partially closed states) "is really impressive"; but f692e
oversplit and some volumes had too little signal. Best combination agreed: f692e discovery (all particles, family
refinement, region boards + arm zooms) + IgM region-stage fixes + gold-standard pruning + polish (minority states frozen).
- Reproduction (run 1): f692e chain from tag f692e-2026-10-18 deployed to routes/toy55_f692e (paths rewritten,
  data symlinked, fboard.sh patched to private project copies), RTAG f692r, job 54020180; watcher f692r_chain.sh
  compares each stage with f692e (f692r_compare.sbatch) and reruns zooms + fullw consolidation.
- Board-level regression test: f692e_repro_test.sbatch (bit-identical, see [[klem-f692e-branch-snapshot]]).
- Pruning test: splitbench8.py (--focus MASK: region mask resampled to the crop grid multiplies every reconstruction,
  null included), prune_prep.py (per-half candidates from level1/states_clusters.txt via simrnasep692_{A,B}/meta_rows,
  FLOOR 500 per half), prune_one.sh (split test only: --pair identity --xem 0 --softem 0 --ctfabs 1, poses =
  poses_famrefined; arm_readout on the board's states_dir), prune_post.py (merged groups, range of motion kept =
  max masked distance among survivors / among all, closed-arm >= 0.8 sub-states kept as own state), prune_all.sbatch
  (LIST, NCHUNK, CHUNK); 75 zooms + 12 region boards, outputs rnasep_klem/f692e_prune/<fam>_<board>/, summary.txt.
  Smoke (fam1 zoomregionmotion_1_s05): 4/4 supported (pair 0.66-0.82 vs null max 0.38 inside a 1.3 % mask), range 1.00.
  CAVEAT: f692e labels came from both halves -> pre-test bias toward support; merges are conclusive, survivals an upper bound.
- PRUNING RESULT (87 boards, 10:56): 382 sub-states, 46 (12 %) below 500 per half, 336 tested -> 332 supported (4 merges),
  range kept 1.00, closed-arm >= 0.8 sub-states in 3 fam1 zooms all kept. As feared: labels from both halves pass almost
  everything; only the floor prunes. -> per-half discovery needed.
- NORM RUNS (2026-10-20 ~11:05, scripts routes/toy55_f692e_norm = the tag set with -p norm, SUBSETFILE half option):
  f692hA / f692hB = chain9 LOCALDEV on one half each (SUBSETFILE simrnasep692_{A,B}/meta_rows.txt, CTFABS 1 PSCALE 1,
  ml_reg=no, NFREE_SUB 1), jobs 54021261/2; f692s55d = the 700k 55-conformer synthetic through chain9 DISTRIBUTED on
  norm (NPARTS 20, NP 20; meta/src from prep_synth55_f692.sbatch in rnasep_klem/f692s55; zooms need CTFABS=0 PSCALE=0),
  job 54021330. Watcher norm_runs_chain2.sh: zooms -> fullw consolidation (halves) / readout55 (synthetic). Next: match
  the halves' consolidated states (reproduced = mutual best match of deviation maps above a null).
- 12:00 relaunch for speed (user: "make the per-half runs distributed ... get things to land faster"): f692hA3 / f692hB3
  distributed (NPARTS 8, NP 8, PARMAX 3; routes/toy55_f692e_unl, which landed on norm anyway), launcher prep_halves_unl.sbatch,
  chains 54024214/5; zooms now chunked (zoom_leaves.sbatch NCHUNK/CHUNK, 4 jobs per family, ZOOM_DONE touched by the
  watcher); watchers f692r_chain2.sh (run 1) and norm_unl_chain.sh (synthetic + halves -> match_halves.sbatch vs f692e).
  Run 1 reproduces Oct 18 closely so far: fam1 region_1 11 states with near-identical sizes, motion_1 7 -> 8 states.
- SYNTHETIC 55 @ 700k through the f692e route (f692s55d, distributed on norm, chain 12:44 = 1 h 43): 24 leaf state maps
  -> 11/55 conformers best-matched (mean corr 0.71; states match whole clades, e.g. 0.71-0.75 to 8/9/10/15/17), 13 zoom
  sub-states -> 9/55 (adds 16, 26) -> 13/55 union; NO closed state: arm readout over all 37 maps max closed-arm 0.59,
  none >= 0.6 (arm_synth55.sbatch). The synthetic closed population (9.3 %) is spread over 7 closed conformers (~1.3 %
  each, different elsewhere) - unlike real data where f692e found a pure closed zoom (0.94). Zoom chunks: 4-6 min each.
Related: [[f692e-validation-cross-half-2026-10-19]], [[klem-polish-loop-refine3d-states-2026-10-20]], [[igm-like-synthetic-armsigm1-2026-10-19]].
