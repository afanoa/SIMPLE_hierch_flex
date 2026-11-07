# Archived refs (2026-11-07)

These refs were on the fork for a few hours on 2026-11-07 and now live only in the git bundle
`/mnt/beegfs/elmlund/afan/mac_archive_20261107/fork_archive_20261107.bundle` on batch2 (114 MB, verified). The bundle is incremental:
restore into a clone that has `origin/master` and `hetero-tree`. Each ref is stored as `refs/remotes/fork/<ref>`:
`git bundle list-heads <bundle>` lists them, `git fetch <bundle> refs/remotes/fork/mac/flex2D:refs/heads/flex2D` restores one.
"Only here" counts commits that were on no remote of the Mac (fork, origin, batch2, box, bw). Uncommitted tracked edits of each Mac
worktree are `git stash create` commits on `mac/dirty/<worktree>` (parent 1 = the branch tip, parent 2 = the index); apply one onto
its branch with `git checkout <branch> && git stash apply <sha>`. Untracked files are in the batch2 tarball (`mac_untracked_list.txt`).
`mac/claude/confident-fermi-0e1cc4` is not in the bundle: its tip is on origin/master.

## Other refs in the bundle

| Ref | Tip | Content |
| --- | --- | --- |
| `wip/state-eval` | `eb58a5455` | `refine=state_eval` and `SIMPLE_STATE_PRIOR` (guide 5.2), built at `hetero_tree_20261028/state_eval` |
| `wip/prob-state-marg` | `8d427a32a` | `prob_state_marg=yes/no` and the `write_states` program; parked, no gain |
| `wip/tree-ports` | `89604e8ac` | the ports, identical to `hetero-tree-wip`'s code |
| `wip/tree-ports-charts` | `480bc9f2f` | `chart_coords=yes` (composed chart coordinates); parked: seams |
| `wip/prior-parentpose` | `df5a3d266` | the pose-anchor session's scripts and verdict |
| `handoff-2026-11-07` | `a8e9070d8` | the first handoff branch: memory, notes, batch2 scripts and 154 MB of the Mac's untracked notes and code |
| `hetero-tree` | `1d474d6e0` | hetero-tree's earlier tip: `e1d5fdd4c` plus the guide's eight doc commits |

## Branches

| Ref | Tip | Not on fork | Only here | Last commit |
| --- | --- | --- | --- | --- |
| `mac/band-exp` | `8588389c7` | 5 | 0 | 2026-10-09 EXPERIMENT (viewing): SIMPLE_FLEXCLS_SWEEP=yes writes cls_expansion_s |
| `mac/batched-insertion` | `bc7889373` | 1 | 1 | 2026-08-11 flex_pca: batch the state-reconstruction plane insertion |
| `mac/claude/amazing-tesla-c2b47f` | `e8da598cd` | 3 | 3 | 2026-11-02 hetero_tree pose-conditioned re-test prototype: a split re-tested on  |
| `mac/claude/confident-fermi-0e1cc4` | `d592c90e9` | 5 | 0 | 2026-11-06 Clean up refine2D/refine3D execution ownership and batching |
| `mac/claude/flex-sort-components` | `bc3a21696` | 2 | 2 | 2026-10-25 flex_pca: each component's Gamma and the mixture prior follow the M-s |
| `mac/claude/objective-meitner-e64b05` | `dc72c27fc` | 1 | 1 | 2026-11-02 hetero_tree: tree_ncores, a core budget for the node steps (budget /  |
| `mac/claude/suspicious-dewdney-bbfa4e` | `998aee001` | 1 | 1 | 2026-11-02 hetero_tree: the per-particle information levers at the arm site prot |
| `mac/claude/tree-mask-discovery-stop-7503a4` | `9b26a8076` | 4 | 4 | 2026-11-02 hetero_tree mask / discovery / stop note and the variance-map prototy |
| `mac/claude/tree-micelle-hetero-73833f` | `8aea8702b` | 1 | 1 | 2026-11-02 hetero_tree: the mic_het fixture, PfCRT in its micelle with large dis |
| `mac/claude/tree-nu-pooled` | `2a4913b42` | 1 | 1 | 2026-11-02 hetero_tree: non-uniform tree-pooled leaf maps prototyped and rejecte |
| `mac/cls-expansion-clean` | `22891675d` | 2 | 0 | 2026-10-09 cls_expansion: conditional sub-class averages with even/odd halves |
| `mac/cls-expansion-mask` | `498ff04e7` | 1 | 0 | 2026-10-14 cls_expansion: soft solvent prior on the class model's basis |
| `mac/cls-expansion-render` | `b66b64f5f` | 36 | 4 | 2026-10-08 cls_expansion: dead code out - the cross-fit (flex_cls_fit_crossed, e |
| `mac/cls-expansion-slim` | `e39834e91` | 2 | 1 | 2026-10-07 cls_expansion note: the ablation table complete |
| `mac/cls-expansion-solvent` | `22891675d` | 2 | 0 | 2026-10-09 cls_expansion: conditional sub-class averages with even/odd halves |
| `mac/coupled-2d3d-rnasep` | `9eaef62df` | 1 | 0 | 2026-10-13 coupled 2D-3D: reconstruction of the coupled-2d3d tree (RAG handoff,  |
| `mac/flex-band-sort` | `90c90a30b` | 24 | 0 | 2026-09-20 flex_pca: a varimax-rotated resume never adopts the cache's deconvolv |
| `mac/flex-cls-split` | `3f687e947` | 87 | 0 | 2026-10-05 TEMPORARY: ablation file switch (nonuis / nocross / noresid) for the  |
| `mac/flex-delivery` | `a94d6c790` | 70 | 0 | 2026-09-12 flex_pca: delete the band-pass fit machinery (SIMPLE_COV_HP and the s |
| `mac/flex-focus-2d` | `cf496b4f6` | 38 | 1 | 2026-09-10 flex_pca: snapshot of the flex-pca-resident working tree (focus mask  |
| `mac/flex-nuis` | `e8ac746c7` | 40 | 0 | 2026-09-29 flex_pca: fixed nuisance block (SIMPLE_COV_NUISANCE=<dir>): generated |
| `mac/flex-pca-merge` | `275844acb` | 12 | 0 | 2026-09-12 wip: worker consumes pindfile |
| `mac/flex-pca-pcg` | `3e80b4265` | 95 | 0 | 2026-09-16 flex_pca: in-engine latent figures and a one-file embedding cache who |
| `mac/flex-pca-reproducibility` | `689577a56` | 17 | 14 | 2026-08-31 refactoring flex_pca to split monolith |
| `mac/flex-pca-resident` | `0fcd70404` | 37 | 0 | 2026-09-09 flex_pca: deflate the consensus dilation (breathing) alongside the en |
| `mac/flex-pca-unified` | `4c7f4bb7a` | 100 | 0 | 2026-09-18 flex_pca master: part scripts request the per-worker nthr, not the bo |
| `mac/flex-polar-em-nuisance` | `9e81476c7` | 14 | 11 | 2026-08-28 Merge remote-tracking branch 'box/flex-pca-em-perf' into flex-polar-e |
| `mac/flex-region-pca` | `1fe238783` | 38 | 0 | 2026-09-20 doc: per-region axes from the 3-column combined embedding |
| `mac/flex-region-pose` | `90c90a30b` | 24 | 0 | 2026-09-20 flex_pca: a varimax-rotated resume never adopts the cache's deconvolv |
| `mac/flex-softblock` | `778458575` | 39 | 0 | 2026-09-27 flex_pca: soft blocks (SIMPLE_COV_SOFT_BLOCKS=<lambda_rel>): a [0,1]  |
| `mac/flex-softem` | `a6b37a2b7` | 9 | 9 | 2026-10-05 softem notes: real-data A/B (10180, RNase P S4), metric calibration,  |
| `mac/flex2D` | `d1a20eff6` | 160 | 0 | 2026-10-03 flex2D: with FLEX2D_CLS_MARG the class bases are deflated of their te |
| `mac/h0901-int-local` | `4ea25d7c8` | 32 | 32 | 2026-09-02 flex: hard-balanced Sinkhorn states inside continuum macro-clusters |
| `mac/state-refs-shrink` | `d1bdc4f96` | 2 | 2 | 2026-10-05 state_refs=consensus: design note and toy/real-data evidence |
| `mac/tmp-port-probe` | `b01c9e88e` | 72 | 0 | 2026-09-16 flex_pca: per-stack sigma2 falls back to the pooled spectrum by itsel |
| `mac/wip-flex-pca-2026-08-11` | `46f9a42ec` | 1 | 1 | 2026-08-11 WIP: flex_pca speedup, auto state discovery, and session experiments |

## Uncommitted work (worktree snapshots)

| Ref | Snapshot | Changed files / branch | Mac worktree |
| --- | --- | --- | --- |
| `mac/dirty/SIMPLE-main` | `09b7c172a` | 20 files on flex2D | `/Users/atifao/SIMPLE` |
| `mac/dirty/SIMPLE_bandsort` | `ba5bc8d7a` | 1 files on flex-band-sort | `/Users/atifao/SIMPLE_bandsort` |
| `mac/dirty/SIMPLE_clsband` | `f65d412a3` | 2 files on band-exp | `/Users/atifao/SIMPLE_clsband` |
| `mac/dirty/SIMPLE_clsclean` | `24487b08e` | 4 files on cls-expansion-solvent | `/Users/atifao/SIMPLE_clsclean` |
| `mac/dirty/SIMPLE_clsmerge` | `4cd79ccc1` | 5 files on HEAD | `/Users/atifao/SIMPLE_clsmerge` |
| `mac/dirty/SIMPLE_flexref` | `5ee1ce195` | 20 files on flex-weights-refine | `/Users/atifao/SIMPLE_flexref` |
| `mac/dirty/SIMPLE_master` | `41840fc95` | 4 files on master | `/Users/atifao/SIMPLE_master` |
| `mac/dirty/SIMPLE_master1013` | `515616f33` | 35 files on HEAD | `/Users/atifao/SIMPLE_master1013` |
| `mac/dirty/SIMPLE_polargpu` | `c61f43351` | 17 files on flex-polar-gpu-2026-08-17 | `/Users/atifao/SIMPLE_polargpu` |
| `mac/dirty/SIMPLE_posemarg` | `42b8a8bf6` | 10 files on flex-pose-marginal | `/Users/atifao/SIMPLE_posemarg` |
| `mac/dirty/SIMPLE_regpose` | `56d023c8a` | 10 files on flex-region-pose | `/Users/atifao/SIMPLE_regpose` |
| `mac/dirty/adaptive-split-merge-b384e9` | `c5db5304f` | 1 files on claude/adaptive-split-merge-b384e9 | `/Users/atifao/SIMPLE/.claude/worktrees/adaptive-split-merge-b384e9` |
| `mac/dirty/codex_worktrees_rnasep-flex-review_SIMPLE` | `19019c4c8` | 1 files on codex/rnasep-flex-review | `/Users/atifao/.codex/worktrees/rnasep-flex-review/SIMPLE` |
| `mac/dirty/epic-rosalind-4bcb43` | `538ad5fe8` | 1 files on claude/epic-rosalind-4bcb43 | `/Users/atifao/SIMPLE/.claude/worktrees/epic-rosalind-4bcb43` |
| `mac/dirty/funny-pike-0046ee` | `0043e3247` | 1 files on claude/tree-dist | `/Users/atifao/SIMPLE/.claude/worktrees/funny-pike-0046ee` |
| `mac/dirty/happy-lichterman-b94a1b` | `7e5367781` | 4 files on landscape-sample | `/Users/atifao/SIMPLE/.claude/worktrees/happy-lichterman-b94a1b` |
| `mac/dirty/intelligent-rhodes-abc9a6` | `2ef724e8f` | 15 files on flex-falloff-prior | `/Users/atifao/SIMPLE/.claude/worktrees/intelligent-rhodes-abc9a6` |
| `mac/dirty/jovial-jennings-821b3d` | `9cdd138aa` | 1 files on claude/tree-data-lp | `/Users/atifao/SIMPLE/.claude/worktrees/jovial-jennings-821b3d` |
| `mac/dirty/objective-meitner-e64b05` | `5abd1b155` | 1 files on claude/objective-meitner-e64b05 | `/Users/atifao/SIMPLE/.claude/worktrees/objective-meitner-e64b05` |
| `mac/dirty/site-states-8c1b66` | `38d5c61a4` | 7 files on claude/site-states | `/Users/atifao/SIMPLE/.claude/worktrees/site-states-8c1b66` |
| `mac/dirty/start-point-sensitivity-cb0d4a` | `74fe5ed7a` | 5 files on claude/start-point-sensitivity-cb0d4a | `/Users/atifao/SIMPLE/.claude/worktrees/start-point-sensitivity-cb0d4a` |
| `mac/dirty/swing-switch-extremes-58027e` | `ee88bd859` | 4 files on claude/swing-switch-extremes-58027e | `/Users/atifao/SIMPLE/.claude/worktrees/swing-switch-extremes-58027e` |
| `mac/dirty/tree-continuum-a56a53` | `cc9d45b0b` | 1 files on claude/tree-continuum-a56a53 | `/Users/atifao/SIMPLE/.claude/worktrees/tree-continuum-a56a53` |
| `mac/dirty/tree-detection-limit-b31e57` | `260c8405d` | 1 files on claude/tree-detection-limit-b31e57 | `/Users/atifao/SIMPLE/.claude/worktrees/tree-detection-limit-b31e57` |
| `mac/dirty/tree-final-routing-f34aad` | `295610c79` | 1 files on claude/tree-final-routing-f34aad | `/Users/atifao/SIMPLE/.claude/worktrees/tree-final-routing-f34aad` |
| `mac/dirty/tree-micelle-hetero-73833f` | `99e5f09fc` | 1 files on claude/tree-micelle-hetero-73833f | `/Users/atifao/SIMPLE/.claude/worktrees/tree-micelle-hetero-73833f` |
| `mac/dirty/tree-micelle-synth-24763a` | `b6ca9863a` | 9 files on claude/tree-micelle-synth-24763a | `/Users/atifao/SIMPLE/.claude/worktrees/tree-micelle-synth-24763a` |
| `mac/dirty/tree-port-optim-de0b3e` | `632164c97` | 8 files on claude/tree-port-optim-de0b3e | `/Users/atifao/SIMPLE/.claude/worktrees/tree-port-optim-de0b3e` |
| `mac/dirty/tree-site-split` | `09b13b997` | 9 files on claude/tree-site-split | `/Users/atifao/SIMPLE/.claude/worktrees/tree-site-split` |
| `mac/dirty/vigilant-cartwright-dcd646` | `c1c92d566` | 16 files on claude/vigilant-cartwright-dcd646 | `/Users/atifao/SIMPLE/.claude/worktrees/vigilant-cartwright-dcd646` |
