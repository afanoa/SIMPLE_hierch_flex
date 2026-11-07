---
name: flexsub-within-state-deformation-no-go-2026-11-07
description: "NO-GO (2026-11-07): subtracting each particle's within-state flex deformation (dz = posterior latent minus its stage-B state's mean, frozen stage-B poses) before recB/stage B/split changes nothing on qtoyflex/qtoy2/qtoynull; per-particle deformation is ~0.05 % of band power; don't retry linear subtraction"
metadata:
  node_type: memory
  type: project
  originSessionId: bf0089df-3038-4376-8273-c1a14df77eba
  modified: 2026-11-07T05:19:43.249Z
---

Tested 2026-11-06/07 on the mm2 gate roots (qtoyflex, qtoy2, qtoynull) with replays from recB (update_project
ptcl_old_root/ptcl_new_root onto corrected stacks, saved recB/stage B/split calls, merge_master2, seed 1, split env). Scripts and
runs: batch2 `math_reworks/latent/flexsub/`; notes `flexsub_within_state_deformation_status.md` in the side-session folder.
- Subtracted component 0.001-0.03 % of the flex-band residual; held-out (checkerboard) gain 0.004 nats/particle on qtoyflex
  (ratio 0.02 of in-sample); within-state R^2 of the arm coordinate s1 on dz 0.02-0.06. The node's own model puts the within-state
  deformation at ~0.05 % of a particle's band power, so the posterior mean is ~96 % shrunk.
- Replays: labels move 0.4-1.1 % (a seed change alone moves 2.9-3.9 %); qtoyflex closed F1 .169 -> .166/.168, holo .619, same arm
  cut; qtoy2 NMI .656 -> .656/.655, K* 5; qtoynull 1 child (split z 4.7 -> 2.1); state-map FSCs move <= 0.0006.

**Why:** in a linear latent model the within-state term is zero-mean around the state mean, so subtraction cannot un-blur a state map,
only remove variance, and that variance is ~0.05 % of the per-particle power. Where dz has signal (qtoy2 states merging conformers,
held-out ratio 0.12) it is discrete sub-state evidence the children need: subtracting it would hurt the recursion.
**How to apply:** do not retry per-particle linear deformation subtraction (or a "perp" variant) before reconstruction; sharpening a
continuum needs a nonlinear warp, which the frozen-pose gate excludes. The uncorrected replay from recB reproduces production bit for
bit, a cheap harness for other per-particle image edits. Related: [[flex-to-refinement-coupling-verdicts-2026-11-06]],
[[flex-frozen-poses-only-no-woodbury-2026-11-06]].
