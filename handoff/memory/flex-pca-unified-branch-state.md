---
name: flex-pca-unified-branch-state
description: 2026-09-16 state of flex-pca-unified (delivery + PCG port): builds at 89e912a53, suite jobs 53228116-19 vs pcgfinal baseline; the hunk-filter port let ~1000 lines of collateral through and had to be redone strictly
metadata:
  type: project
---

`flex-pca-unified` (Mac worktree /Users/atifao/SIMPLE_unified; batch2 staging ref
`flex-pca-unified-in`, build tree `~/SIMPLE_unified_p2`, script `build_p2.sh`) = flex-delivery
a94d6c790 + master d8ae11d15 + the PCG line ported as hunks + figures/one-file cache + defaults.
First `make rc=0` on 2026-09-16 15:46 at 89e912a53; three runtime defects fixed by 59d9769f5
(sigma fallback not in part scripts -> set_worker_key; plain v5 part reader vs band-boxed writer;
paired path returned before the state stage -> l_paired_states). Suite tag `unified3` on
59d9769f5; baseline `pcgfinal` (jobs 53227914-17, `~/SIMPLE_pcgfinal`, tag flex-pca-pcg-final =
3e80b4265) scored ribo 5/4/K4, splice 5/2/K5, trpm4 8/5/K8, cnga1 3/2/K4 -- ribo and splice
UNDER the reference table (6/2/K5, 10/3/K16) on the same recipe: unexplained as of 09-16. Score with
`python3 ~/flex_verify/tools/suite_score.py pcgfinal unified`. Frozen verified tree:
`~/SIMPLE_pcg_opt10` / tag flex-pca-verified-2026-09-16 (leave alone).

**Why:** the regex hunk filter used for the port (keep hunks containing a PCG identifier) let
collateral through wherever pcg-only motion shared a hunk: pcg's whole `run_flex_pca` body and
state delivery in model.f90 (about 1000 lines with no PCG identifier), pcg's per-state delivery
in rec3D, focus-mask and DEFLATE_CTF lines in mstep/pairmerge, polar imports in the em parent,
duplicated bodies in em_fit. model.f90 needs NO PCG code at all.

Two more port gaps found by an iteration-1 log diff against pcgfinal (09-16 night): the per-fit
M-step's DEFLATE_BG defaulted off (now on, =0 opts out, b85a5cf7c) and the envelope window was
missing at two delivery-side sites (realize_hermitian_volume: sphere instead of envelope -> the
calibration Gram 4.0x off; merge.f90 state masking), 3731ffd13. Suite tags: unified3 (no BG),
unified4 (BG, sphere reps), unified5 (both fixes) -- unified5 REPRODUCED pcgfinal's basis
(iteration-1 numbers to 3-4 digits; merged components splice 7/2 vs 5/2, ribo 5/4 vs 5/4;
principal angles identical to 3 decimals). unified6 (e8cd25997, `~/SIMPLE_unified_p3`) adds pcg's
rec3D state delivery (per-state own-kernel gridding, own eo-FSC low-pass) replacing the delivery's
per-shell pooling; cnga1 distributed parts need MEM=240000. unified6 REPRODUCES pcgfinal end to end on splice/ribo/trpm4
(states 9.7-8.9 / 4.2-4.9 / 7.8-8.5 A vs 9.7-8.9 / 4.2-4.9 / 8.0-8.7 A); the both-halves-at-0.7 count
is noise-sensitive on trpm4 (8 vs 4 with identical iterations) -- read meanmin and the per-iteration
agreement too. cnga1 at neigs=10 (`unified6k10`) is the candidate recipe (K=32 = 528 pair kernels
is the whole cnga1 cost).

**How to apply:** the fast port audit is a line diff of the first EM iteration's log against the
reference run: EM CALIBRATION (sig2, gamma0, tr(G)), deflation rank and energy, |b|, |x0|,
-2logL/N. Any factor that is exactly a support-volume ratio means a sphere/envelope mismatch.
 when porting between the delivery and pcg layouts, reset the file to the base
and re-add named PCG items by hand (em parent: two use lines, probe_fit_t members, the kpk/rpk
arguments of the four part-file interfaces; em_iter: four hooks; rec3D: the rec_states_backend
dispatch only). Audit with a per-hunk pcg-line count before building. `ssh batch2` (short alias)
now fails host-key verification in BatchMode; use `batch2.ncifcrf.gov`. See
[[flex-pca-fork-map-2026-09-15]], [[flex-pca-pcg-branch-state]].

**2026-09-17 memory/speed (exact):** band-ball index lists in flex_pcg_t (0eca9f824..8d030b8e1):
cnga1 K=32 worker peak 197 -> 37 GB, worker round 648 -> 319 s; kernels bitwise (test D), pipeline
within its 1e-3 run-to-run band. Reduce unchanged (parts were already index-listed on disk). The
finalize was strided-store bound; slab + blocked transpose. Test equivalence on a 5k-particle
subset (`pair_test.sh` with NPTCL=5000, ~4 min), never the full set -- the user objected to the time.

**Squash (2026-09-17):** `flex-pca-memspeed` = single commit 2cc35d3ac on d8ae11d15 (48 files, no
docs), tree = flex-pca-unified tip minus 4 docs; for the merge to master (5 files conflict: model,
ui_denoise [master moved flex_pca UI to simple_ui_heterogeneity.f90], flex commander, phases, test).
Master's refine3D_states flex=yes path needs min_state_frac + place_states_with_population_floor +
pickup_project_consensus_volume brought across; our branch already writes state labels and
registers 'vol_flex' at native box. cnga1 K=32 end to end: 118 min vs 162, same scorer row.

**MERGED TO MASTER 2026-09-17:** 8b9922b5c "half set EM engine, pcg kernels, solvent mask support,
parameter tidying" (one commit, 48 files, +16,953/-10,139, fast-forward on 03b3a098b). Master's
population floor (min_state_frac) and consensus-map pickup were ported by hand; flex_pca UI now lives
in simple_ui_heterogeneity.f90. Second commit 3a4c7412e (refine3D_states flex=yes: labels only,
rec_states=no, gridding) is on `flex-pca-memspeed` in /Users/atifao/SIMPLE_landing, built (make
rc=0), NOT pushed -- user's call. flex-pca-unified keeps the 98-commit history.
