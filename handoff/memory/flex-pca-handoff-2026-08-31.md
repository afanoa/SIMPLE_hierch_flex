---
name: flex-pca-handoff-2026-08-31
description: "three-stream handoff (half-set EM+rank, GPU plan, common-lines volume-free M-step) and the branch-split doc-location trap"
metadata: 
  node_type: memory
  type: project
  originSessionId: 74eda56d-80ed-4b21-b788-b0776c91d800
  modified: 2026-08-31T19:57:59.844Z
---

Handoff doc: `doc/for_developers/ideas/flex_pca_handoff_2026-08-31.md` (untracked on the
master worktree as of writing). Three streams, each with a proposal + pre-registered Phase 0:

1. **Half-set EM + varimax rank** — `flex_pca_halfset_em_varimax_rank_proposal.md`; see
   [[flex-pca-varimax-rank-criterion]]. Merge-don't-refit final stage; FSC doctrine (3 roles,
   Gilles–Singer RECOVAR SI at `~/Downloads/pnas.2419140122.sapp.pdf`); sampling-aware
   S.11/S.13 ridge; `flex_pca_crossfsc.bin` carries curves + sampling H.
2. **GPU/CPU** — `~/ribo_analysis/GPU_CPU_PLAN.md`, measured on Blackwell: end-to-end 1.75x,
   M-step insert GPU 1.2 s vs CPU 11–30 s (anti-scales past ~11 threads; z-slab rethread in
   worktree is the fix), STATEREC prep ~45 s dominates both arms, fixed tail caps ~2x.
   P0 ship measured / P1 z-slab / P2 staterec prep / P3 window-safe E-step. Do-nots: bank
   tuning, full-pass GPU E-step, GPU states/merge.
3. **Common lines / volume-free M-step** — `flex_pca_volume_free_mstep_proposal.md`
   (IUCrJ 2026, 13, 354–363 basis): fan ring-algebra M-step + per-iteration band-limited
   realization; ring FRC under a shared basis is INVALID (E[z·n] common-mode — half-set EM
   is a precondition); direction quantization not curable by ndir; 3–6x on the rank ladder.

**The doc-location trap:** the Aug 26–28 corpus (47 docs incl. the old halfsets proposal,
staged plan, speedup review, speed_and_quality findings, polar_em proposal, 08-24 handoff;
commits `7a6d9acb0` + `cc78d6249`) lives ONLY on branch `flex-pca-reproducibility` — NOT in
master's tree. Read via `git show flex-pca-reproducibility:doc/...`. Current checkout moved
to master @ `5aee28fe7` (split-monolith lineage, on origin; local master 17 behind
origin/master on 2026-08-31). New proposals are untracked on master.

**How to apply:** start any flex_pca session from the handoff's §6 ordered actions (it_win
fix first); when a cited doc is missing from the tree, check the branch before concluding
it was deleted.
