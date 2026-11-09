---
name: wip-settling-pass-prob-beta-cold
description: deterministic settling of wip refine3D without code: one stage-B iteration with prob_beta=1e9 collapses every keyed draw (pose, in-plane, shift disk) to its mode (n_eff 1.00 vs 2.53 at beta 1)
metadata:
  type: project
---

On the wip tip (918215446 and the nodelevers test build) all three keyed draw families of the 3D search are drawn at
prob_beta x T and the draw routine returns the mode when the temperature is effectively infinite, so appending the stage-B
refine3D call once with `maxits=1 prob_beta=1e9` from the final maps is a posterior-mode pass: the log reports
"POSTERIOR DRAW N_EFF: 1.00 at beta 1e9, 2.53 at beta 1.00" (r1, 252k particles, 2026-11-09). prob_beta is an existing CLI key
(stage A uses 0.3); its other uses are readouts and neff logging. Test harness: node_levers/nl_settle.sbatch (+ nl_settle_summary.sh).

Tested on the six arm-b r1 replays (2026-11-09): the first cold pass moves 58 % of directions, a second 5.2-5.7 % (fixed point);
hit / miss unchanged, settled-label ARI +.03 within a basin (.74-.77) and .26 across; matched-state responsibility correlation .95
within basins either way; the closed child's map gains closed-arm +.016 / open-arm -.03.

**Why:** the wip refinement re-draws 72 % of directions per iteration ([[wip-refine3d-redraws-poses-each-iteration-2026-11-08]]);
the user wants stored poses and labels deterministic, and no new switch was needed.

**How to apply:** for deterministic final poses/labels or for per-particle state evidence at the mode, run the cold pass; do not
rebuild for a settle switch. The test build's writer (SIMPLE_STATE_EVIDENCE=yes) then gives the responsibilities at that mode.
