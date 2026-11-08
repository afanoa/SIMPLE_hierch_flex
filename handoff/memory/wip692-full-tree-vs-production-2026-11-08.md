---
name: wip692-full-tree-vs-production-2026-11-08
description: "USER 2026-11-08 'Yes, run it tonight': one full 692k RNase P tree on hetero-tree-wip defaults (master merge + depth_wiener + noise model default-on) vs production's last 692k run f692_lf_s1 and r17d, scored by production's own scorecard with better/same/worse rules fixed before the run (guide 7.10 thread 11)"
metadata:
  type: project
---

The user asked whether anything developed on wip improves RNase P's closed state over production. Subtree and fixture tests say no
once read honestly, and no wip build had run at 692k. So: one full tree tonight.
- Session wip692 (local_facb8c0b-c228-414a-b0c8-fc41ce41e8c8). Run `$H/integration/f692_wip_s1`, notes `$H/handoff/wip692/STATUS.md`.
- Build: rnasep-misfit's nvon build (wip 0a29cfbee + noise default-on), with the OpenBLAS variant matched to build_lf. Launcher
  and keys are exactly f692_lf.sbatch's.
- Baselines: f692_lf_s1 (5 h 14, 617 core-h, closed leaves r1312_g1_2 / r2112_g1_2 / r1321_g1_2, purest .095, union F1 .137;
  landscape arm CC .800 / .705) and f692_r17d_s1 (5 h 31, 873 core-h, corner r1311_s2, purest closed .094).
- Rules fixed before the run:
  - Gate: at least one leaf the user judges closed.
  - Same: inside the lf / r17d span.
  - Better: closed kept, plus purest closed leaf >= .11, best closed-leaf F1 >= .14 or union F1 >= .16, with the landscape no
    worse and wall <= 7 h.
  - Worse: fewer accepted closed leaves, purity < .08, F1 < .10, a shrinking landscape, or wall > 7 h.
- Caveats: one seed, own-data maps. A clear result needs a second seed and a production-tip twin (merge_master2) to attribute it.
**Why:** the first apples-to-apples full-scale test of the past days' code against production.
**How to apply:** judge the result by these pre-set rules, closed reference rows first; never move the bar after seeing it.
Related: [[integration-692k-lf-launch-2026-11-02]], [[integration-lf-closed-three-leaves-2026-11-02]],
[[rnasep-closed-state-must-be-retained]], [[rnasep-misfit-noise-spectrum-verdict-2026-11-08]].
