---
name: tree-continuum-fixture-and-verdict-2026-10-29
description: "hetero_tree on a 1-D continuum (ctpk dwell / ctun uniform swing toys): baseline stops after 1 cut; SIMPLE_TREE_INHERIT (default off) reaches the closed end; occupancy unfixable via counts because poses push leaf maps outward"
metadata:
  node_type: memory
  type: project
  originSessionId: f6f31e73-5d9d-4609-b06e-14cf14e50c04
  modified: 2026-10-30T00:38:47.056Z
---

2026-10-29, session tree-continuum. Work dir: batch2 hetero_tree_20261028/continuum/ (scripts/, sim_ctpk, sim_ctun, runs/,
diag/, RUN_ENV.txt). Shared-notes write-up: continuum_2026_10_29.md + inherit_on_r17d.diff (md5 01764c0f).

- Fixture: scripts/run_ctoy_sim.sbatch = qtoynull recipe + rigid swing of the closed conformer's domain (make_swing_maps.py, 30 A tip,
  31 x 1 A steps). ctpk = dwell at t 0.65 (564/step at the closed end, 5,244 peak); ctun uniform. Truth t = closed_arm column.
  Scorer scripts/ctoy_score.py (member t + NNLS map content; the best-step CC is biased for mixtures).
- Per-particle position at qtoy SNR: sd 0.33 of the span (best third 0.27). Evidence tails cannot isolate the ends.
- Baseline recipe17d: root cuts once at the mass median; children have no pole -> fflex fallback seeds an unrelated axis -> FAIL.
  2-3 leaves, no map in t 0.15-0.7, far end 0.76-0.78.
- Leaf maps sit 0.1-0.3 further out than their members: pose error grows with distance from the map's own end (true-pose
  reconstruction = members). So leaf counts are not lifetimes; soft mixture weights over the maps do not recover occupancy either.
- SIMPLE_TREE_INHERIT=yes (build_inh): a pole-less node below the root splits at the median of its parent's stage-B evidence
  contrast. 4 leaves; closed-end maps 0.00-0.08 in 4/4 runs; far end 0.81-0.88; ctun reproduces across seeds, ctpk interior doesn't.
  qtoy2: changes only pole-less r1 / r4 (NMI 0.686 -> 0.685, closed unchanged, holo over-split).

**Why:** user goal = leaves that produce the extremes and the correct sampling of a continuum.
**How to apply:** default stays off. The fflex fallback seeds 71/84 nodes of the 692k tree, so adoption needs the RNase P closed
gates + OLE first. Related: [[locality-guard-skips-em-group-splits-2026-10-29]], [[import-particles-resets-state-and-eo]],
[[pose-state-coupling-toy-verdict]], [[split-test-accepts-quality-split-on-null-2026-10-27]].
