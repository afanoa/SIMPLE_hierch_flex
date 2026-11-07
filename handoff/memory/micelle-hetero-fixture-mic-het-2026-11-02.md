---
name: micelle-hetero-fixture-mic-het-2026-11-02
description: "mic_het (hetero_tree_20261028/fixtures/pfcrt_micelle/mic_het): PfCRT micelle fixture with LARGE heterogeneity (TM half rocked 35 deg, Fv swung +/-40 deg, Fv absent; 40/25/20/12.5/2.6 %), matched to micnull. VERDICT: the tree finds 3 clean + 1 mixed state maps (NMI .14, memberships .45-.70 pure, rare lost); the pose flag dflip>=.04&z>=5 FIRES on the real root split (Fv states re-align worse: .40 vs .27 turned, dflip .13) -> not selective, report only, never gate; fixposes after the pose gate recovers 4/5 states (NMI .33) but only because the inputs are 1.4 deg from truth"
metadata:
  type: project
---

2026-11-02, session tree-micelle-hetero (note doc/implementation_notes/tree_micelle_hetero_fixture_2026_11.md on branch
claude/tree-micelle-hetero-73833f; scripts scripts/hetero_tree_micelle_fixture/ and W/scripts on batch2, W = hetero_tree_20261028/fixtures/pfcrt_micelle).
- Generator `micelle_fixture.py maps --spec ref,rock35,fvswing40,nofv,fvswing-40` (ops rock / fvswing / nofv; `--states` is the old
  rock-only form), `oritab --micbins micnull/truth.txt` copies the per-particle micelle bins; `make_michet.sbatch` runs the whole
  fixture in one 16-core job (26 min) and copies micnull's measure.json / measure_cn.json so the noise is identical.
- States: 22,818 ref / 14,052 rock35 (1,506 TM atoms, 10.4 A mean, 22 A max) / 11,312 fvswing40 (1,780 Fv atoms, 13.4 A mean, 29 A max) /
  7,151 nofv / 1,475 fvswing-40 (2.6 %, below the known root detection limit n/sqrt(N) >= 16: here 6.2). Input src_in_fp.simple at the
  shared a3_s0.75 poses (1.41 deg / 0.92 px), no refinement. Truth: mic_het/truth_full.txt. Realism 0.144/0.097/0.061/0.052.
- TRAP: the micelle is carved around the UNION of the state masks (21 % more cavity in the slab than ref alone), so the TM rock is
  protein-vs-vacuum contrast, over-stated; the Fv states live above the belt and are not.
- Trees: trees/michet_head_s{1,2}, micnull_head_md1_s1 (build_head = c7ac3b0bd), michet_gate_nomask_s1 (POSEGATE 0.165 -> fixposes);
  keys = gate keys with minpop 1000, maxdepth 1. Scorer `micelle_score.py --maps` (pair V + z, NNLS map composition), `michet_tables.py`.

- RESULTS (both seeds identical to 160 particles, 15 core-h per tree): root split Fv-in-place vs Fv-displaced (V 0.47) is REAL and a
  pose partition at once (true flips 0.27 vs 0.40; dflip 0.13 z 38-40 -> FLAG); below it ref vs rock35 (V 0.48-0.51, dflip 0.038 /
  0.003, unflagged: seed 1 within 0.003 of the line) and fvswing40 vs a fvswing40/nofv mixture (V 0.16-0.19). Leaf maps: ref .91,
  rock35 .73-.78, fvswing40 .91-.93, nofv only .39 in a mixed map; rare state in no map (n/sqrt(N) 6.2). micnull: one leaf (root
  pose candidate dflip 0.071 flagged, guard merged it as diffuse, L/L_ref 0.24 vs 1.36 for the real split).
- fixposes response (recipe17d_micsup_fpg, POSEGATE 0.165 fired at 32.3 %): NMI 0.33, rock35 / fvswing40 / nofv found (purity
  .83/.63/.55), ref .79 (oversplit, V 0.17), rare 52 % recalled inside the nofv leaf; harmless to real states ONLY because the input
  poses are 1.4 deg from truth (real PfCRT inputs are at chance between refinements).

**How to apply:** use mic_het (with micnull as its null) for any question about real states in a pose-unstable, micelle-dominated
dataset. Never gate on dflip: a real state whose particles re-align worse against the node's consensus produces the artefact
signature; read the 692k r4 (dflip 0.106) with that caveat. The locality of the difference (guard) told the two apart here, without a mask. Related: [[pfcrt-micelle-fixture-tree-pose-scatter-2026-10-29]],
[[tree-pose-artefact-test-2026-11-02]], [[synthetic-toy-traps-noise-poses-2026-10-29]], [[rare-state-detection-limit-2026-10-29]].
