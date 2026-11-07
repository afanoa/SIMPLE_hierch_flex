---
name: klem-fortran-home-state-weight-set
description: "Where the kernel-latent EM goes in Fortran: master's state_weight_set (per-state fractional particle weights, KERNEL kind, parent-state genealogy via publish_work_state, consumed by refine3D/reconstruct3D/rec3D service/eul_prob_tab) + the soft-EM E-step on branch flex-softem (simple_softem_weights.f90, 85185be4c); found 2026-10-15 at the user's hint"
metadata:
  type: project
---

User 2026-10-15: "There might be implementations in other branches/chats you may find useful" and "maybe there are
genetic approaches ... these volumes' parents are these states".

What exists (checked with git grep across local branches, 2026-10-15 18:40):
- `src/main/project/simple_state_weight_set.f90` + `src/fileio/simple_state_weights_file.f90` are ON MASTER (introduced in
  19f458193 "forgotten files"; 23 files on origin/master use it; NOT on the older flex2D branch). A transactional set of
  per-state particle weights (one file per state + manifest; kinds PARTITION = rows sum to 1, KERNEL = fractional), with
  producer / generation / layout digest and a PARENT: `publish_work_state(parent, state, parent_rows, ...)` publishes a
  work project's single-state set derived from one state of a parent set (refine3D_auto state=X builds it) and records
  parent_generation / parent_state = the genealogy of states. Consumers: simple_commanders_rec(_distr), refine3D
  commander and strategy, rec3D service / pcg strategy, eul_prob_tab, matcher_3Drec, flex_pca state service / gateway /
  application, sp_project_out (imgkind state_weights), final_rec.
- Soft-EM E-step in Fortran: branch `flex-softem` (worktree .claude/worktrees/suspicious-sammet-9e8bee, commit 85185be4c
  "softem: per-iteration soft state weights in refine3D (SIMPLE_SOFTEM) and consensus+deviation references
  (SIMPLE_SOFTEM_REFS=sx)"): `simple_softem_weights.f90` = cross-half (pftc eo flip) image-likelihood posteriors over
  states after each refine3D batch search, logsumexp over in-plane x directions, fractions fixed; batch2 tree
  ~/SIMPLE_softem. See [[flex3d-soft-em-toy-results-2026-10-05]] for its findings (references are the lever; drain traps).

**How klem maps onto it:** cells = states (49 + junk); responsibilities = KERNEL weights; the grid kernel = a smoothing
of the weight columns across neighbouring states before reconstruction (the one new piece); M-step = reconstruct3D with
the smoothed weights and consensus+deviation references (SIMPLE_SOFTEM_REFS=sx already does consensus + deviation);
E-step = softem_weights with the cross-half pftc scoring; schedule = lp 20/14/10 with the kernel width annealed; the
hierarchical zoom = publish_work_state children from a parent state's rows (the "genetic" expansion the user means),
with a cross-half reproducibility test deciding which children survive (the 3D analogue of cls_expansion's bisection).
**How to apply:** port after the Python prototype settles (hierarchical zoom validated on toy10 first); build on a
master-based branch (state_weight_set), cherry-pick the softem E-step. Related: [[klem-kernel-latent-em-2026-10-15]],
[[stored-weights-pipeline-toys-then-not]], [[prototype-in-python-before-fortran]].
