---
name: flex-soft-blocks-branch
description: "2026-09-28: branch flex-softblock (Mac worktree /Users/atifao/SIMPLE_softblock, batch2 ~/SIMPLE_softblock, staging ref flex-softblock-stage) = flex-region-pca + SIMPLE_COV_SOFT_BLOCKS=<lambda_rel>: a [0,1] order weight as pcg_mskfile -> protein + nuisance blocks on the sphere, real-space M-step penalties lambda(1-w)/lambda w; toy says lambda_rel 0.3-1"
metadata:
  type: project
---

Commit 778458575 (one file, simple_flex_pca_pcg.f90, +160): module state l_flex_soft / flex_pen_img / flex_seed_img;
flex_pcg_t%pen_blk/pen_mean/lam_pen/pen_rel/l_pen + set_block_penalty; the penalty is added per column in
apply_operator, its block mean folded into apply_precond's diagonal, lam_pen = pen_rel x lam_base derived in
prep_floor; soft_blocks_from_weight (contained in flex_pcg_support_volume, covariance-box store call only) makes
flex_nblk=2 with both block windows = the soft mskdiam sphere, flex_window_apply_rec_region uses the seed windows
(w, 1-w) x sphere (the data-free initialiser is its only caller). Syntax-checked locally against
/Users/atifao/SIMPLE/build/modules (baseline and edited clean); built on batch2 via sbatch ~/build_softblock.sh
(x6342, 4 min, rc=0). Updates: push to flex-softblock-stage, `git -C ~/SIMPLE_softblock reset --hard
flex-softblock-stage`, sbatch the build (never while jobs run from that tree).

Weight input: `agnostic_mask.py` writes `agnostic/<ds>_protein_weight.mrc` = logistic((lcc-thr)/sd_disordered)
inside the belt slab, 1 elsewhere in the automask, 0 outside, blurred 1.5 A (cnga1 mean 0.90 protein / 0.14 belt;
PfCRT 0.87 / 0.16). Toy `toy_micelle/lin_soft.py` (belt touching the protein, 24.6 % of its energy on protein
voxels): leak/recovery A,B,C -- hard protein-only 0.52 / .55 .52 .38; hard two-block 0.16 / .86 .62 .46;
soft 0.3 0.22 / .73 .75 .84; soft 1 0.30 / .70 .66 .82; soft 3 0.52; soft 10 0.71 / .11 .41 .19 (w in the belt is
~0.14 not 0, so strong lambda penalises both blocks there). Soft wins only on the boundary-straddling mode C.
Real-data arms `flex_region/s_hard`, `s_soft0.3`, `s_soft1` (nparts=1, 32 cpu, cn001); score with
`regress_soft.py s_hard s_soft0.3 s_soft1`. Related: [[cnga1-micelle-leak-is-negligible]],
[[agnostic-order-protein-mask]], [[flex-region-pca-branch-state]].

**REAL-DATA RESULTS 2026-09-28.** CNGA1 (nparts=1, same binary): hard order mask leak R^2 0.9 %, 2/10 reproducible
(min|cos| mean 0.61); soft 0.3 leak 3.5 %, protein cols 60-85 % on protein, 5/10; soft 1 leak 1.7 %, 73-91 %, 7/10
(0.75), nuisance block ~ R2 belt R^2 0.45; the soft protein block has LESS boundary-zone variance (0.16) than hard
(0.25) -- the toy's boundary benefit did not appear. PfCRT (16.8k, box 100, neigs 6): hard 99.3 % of protein-block
variance on protein, 0/6 reproducible (homogeneous sample, as expected), leak ~0.1-0.4 %; soft arms spill 15-31 % of
protein-column energy OUTSIDE the automask into solvent and their nuisance blocks are 62-82 % solvent -- on a
signal-free set the soft penalty lets noise roam. Next fix if pursued: union = automask (not the mskdiam sphere, which
also came out at 0.50 of the box on cnga1 via mask3D_soft's edge) so neither block can wander into bulk solvent.
Docking (models/dock.py, 7LFT/6UKJ): agnostic mask covers 99.8 / 99.7 % of model voxels, belt holds 0.0 / 0.2 %;
our cnga1 map is the MIRROR hand of 7LFT (CC 0.77 vs 0.62); the old PfCRT 100 kDa mask misses 11 % of chain A.
PfCRT hard components split energy transporter 18-32 % / Fab 8-20 % / protein shell ~50-68 %: to isolate the helices
the Fab must be excluded too (order cannot tell Fab from transporter; both are ordered).
