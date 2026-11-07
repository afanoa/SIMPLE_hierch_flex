# mic_het: micelle fixture with large heterogeneity (session tree-micelle-hetero, worktree tree-micelle-hetero-73833f, branch claude/tree-micelle-hetero-73833f off claude/tree-dist c7ac3b0bd)

Work dir on batch2: /mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/pfcrt_micelle (W): mic_het/ (dataset), trees/michet_*, trees/micnull_head_md1_s1,
trees/pose/ (pose statistics), scripts/ (micelle_fixture.py extended, micelle_fixture_v1.py = the 2026-10-29 version, make_michet.sbatch,
submit_michet.sh, score_michet.sbatch, pose_michet.sbatch, pose_split.py (+ truth columns), michet_tables.py). Repo copies in
scripts/hetero_tree_micelle_fixture/; note doc/implementation_notes/tree_micelle_hetero_fixture_2026_11.md (uncommitted: the maintainer commits).

## 2026-11-02 ~12:30 (batch2 time)
- Fixture mic_het built (job 54857068, 26 min x 16 cores): micnull's micelle bins, quality draws, noise realisation and band powers copied
  (matched pair; only the protein differs). States ref / rock35 (TM half 35 deg, 10.4 A mean) / fvswing40 (Fv 40 deg sideways, 13.4 A
  mean) / nofv (Fv absent) / fvswing-40 (rare) at 40 / 25 / 20 / 12.5 / 2.5 % (22,818 / 14,052 / 11,312 / 7,151 / 1,475). Input
  src_in_fp.simple at the shared 1.41 deg / 0.92 px perturbed truth poses. Realism 0.144 / 0.097 / 0.061 / 0.052 vs real 0.151 / 0.099 /
  0.069 / 0.069 (the states lower it; the noise is micnull's).
- Trees launched (jobs 54859817-54859820): michet_head_s1 / s2 and micnull_head_md1_s1 on math_reworks/build_head (= c7ac3b0bd, gate 15/15),
  michet_gate_nomask_s1 on recipe17d_micsup_fpg with SIMPLE_TREE_POSEGATE=0.165 (the fixposes response). Keys = the micelle gate keys
  (nparts 5, mskdiam 200, site layer + states, kway_maxdepth 2, CONTRAST=maps, KEEPMAPS) with minpop 1000 (the gate's 3000 = 5.3 % could
  never hold the 2.5 % state), maxdepth 1, nodepar 4. Scoring + pose statistics follow automatically when they finish.

## 2026-11-02 ~13:10 (batch2 time): RESULTS (note sections 3-4 written; commit pending the last control)
- Both seeds identical (160 ptcls): root = Fv-in-place vs Fv-displaced (V 0.47, REAL) and at once a pose partition (true flips 0.27 vs
  0.40): dflip 0.13 z 38-40 -> the proposed flag FIRES on a real split at both seeds. Depth 1: ref vs rock35 (V 0.48-0.51) unflagged
  (dflip 0.038 / 0.003, seed 1 within 0.003 of 0.04); fvswing40 vs fvswing40/nofv mixture (V 0.16-0.19) unflagged. micnull: one leaf
  (root candidate dflip 0.071 flagged, guard merged as diffuse 0.24; the real split is localized 1.36). NMI 0.14; leaf maps ref .91,
  rock35 .73-.78, fvswing40 .91-.93, nofv .39 (mixed); rare 2.6 % state lost (limit).
- fixposes response (fpg build, gate fired at 32.3 %): NMI 0.33, 4/5 states found (purity .55-.83), rare 52 % inside the nofv leaf,
  one ref oversplit; harmless to real states because the inputs are 1.4 deg from truth. Control michet_nomask_s1 (free poses, same old
  build) running to separate build from poses.
- VERDICT: the flag is not selective (real states that re-align worse against the node's consensus look like pose artefacts): report
  both numbers, never gate; 692k r4's 0.106 does not by itself make it an artefact. Cost ~55 core-h so far.
