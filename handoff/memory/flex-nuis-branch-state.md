---
name: flex-nuis-branch-state
description: "branch flex-nuis (e8ac746c7, off flex-softblock): SIMPLE_COV_NUISANCE=<dir> fixed nuisance block marginalised by a Schur fold in the E-step + embedding, Psi EM-updated; batch2 worktree ~/SIMPLE_nuis (own build); arms n_*/q_* on cnga1/PfCRT"
metadata:
  node_type: memory
  type: project
  originSessionId: 2590d10f-7aff-4aff-b2cc-d9d85f962ad3
  modified: 2026-09-29T14:31:13.726Z
---

2026-09-29. Mac worktree /Users/atifao/SIMPLE_softblock now on branch flex-nuis (commit e8ac746c7 on top of the soft-block
commit 778458575); pushed to batch2 as flex-nuis-stage, worktree ~/SIMPLE_nuis (branch flex-nuis-wt), built by
`sbatch ~/build_nuis.sh` (x6342, ~4 min, done-file make_nuis.done). ~/SIMPLE_softblock stays the soft-block binary.

Design (all in src/main/flex, no new files): module state nuis_recs/nuis_psi0/nuis_psi_run/nnuis/l_nuis in the em parent;
probe_fit_t%nu_* per-fit Psi + scratch; estep: nuis_requested/nuis_setup/nuis_fit_begin/nuis_fit_reduce/nuis_plane_stats/
nuis_schur_fold + private nuis_estep_fold (end of the Cartesian former) and nuis_estep_post (after the solve: E[w],
E[ww'] -> Psi accumulator, T E[w] added to mean_fpl so the residual is y - aP(mu + T E[w])). S = G_ww + (sig2/a^2)Psi^-1,
G/b/c -> Schur-folded, then the UNCHANGED probe_solve_ecm/probe_solve_mix (exact at fixed contrast; nml_plain is 0).
Plug-in M-step (Cov(w,z) not inserted). Psi = mean E[ww'|y] each iteration (SIMPLE_COV_NUISANCE_EM=0 keeps the file
prior); paired fits' Psi averaged -> polish (fit id 0) -> embedding (folds full + half solves at a=1 / per-half ah;
writes flex_pca_nuisance_latents.bin: int32 n, nn, pinds; real32 w(n,nn) column-major). The block forces the Cartesian
former (no polar bank slots); SIMPLE_COV_ESTEP_CART=1 = same former without the block (baseline). nparts=1 and CPU only
(nuis_setup throws otherwise). project_fplanes_mean_basis gained optional extra_recs/extra_fpls.
Template dir format = a finished run's delivery: flex_pca_merged_pcNNN.mrc at box_crop + flex_pca_probe_merged.txt
(count, sig2, prior variances). Generator: flex_region/nuis_model2.py <cnga1|pfcrt> (nuis_cnga1, pfcrt/flex_soft/nuis_pfcrt).
Arms: cnga1 flex_region/n_{hard,auto}_{cart,nuis} (launch_nuis.sh), PfCRT pfcrt/flex_soft/q_{hard,auto}_{cart,nuis}
(launch_pfnuis.sh); driver scratchpad/nuis_driver.sh keeps <= 3 jobs; scorer flex_region/score_nuis.py <ds> <ref> <arms>.
Related: [[micelle-estep-marginalization-test]], [[flex-nuisance-block-code-map]], [[flex-soft-blocks-branch]].

**VERDICT 2026-09-29 (full EM, Cartesian former in every arm):** the marginalised micelle block is a clean negative.
Automask + block: PfCRT leak 13.2 % -> 8.9 %, belt still 53 % of component variance (compose test on cnga1: 14 % -> 4.7 %).
Order mask + block vs order mask alone: PfCRT leak 0.08 vs ~0 %, 1/6 vs 0/6 reproducible (mean 0.551 vs 0.545);
cnga1 1/10 vs 3/10 reproducible (mean 0.548 vs 0.641; polar s_hard 2/10, 0.614). The block's E[w] holds 29-31 % of the
belt reference but buys nothing on top of the order mask. User concluded "the order mask is the way to go"; automask arms
cancelled. Principled framing given: heteroscedastic factor analysis -- with protein/belt variability independent (~1 %
shared) its ML factors vanish on the belt, i.e. the order mask IS the FA solution here.
