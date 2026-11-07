# hetero_tree: testing and handoff guide (2026-11-07)

Everything needed to pick up the hetero_tree work (SIMPLE's 3D heterogeneity tree, its coupling with flex_pca, the open experiments)
without the Mac it was developed on. Written for the user and for a fresh Claude Code session. Companion records:
[latent-refinement verdicts](hetero_tree_latent_integration_2026_11_06.md),
[flex -> refinement tables](hetero_tree_flex_to_refinement_2026_11_06.md),
[M-estimation and ML formulations](hetero_tree_mest_formulations_2026_11_06.md).

## 0. Start here

1. `git clone git@github.com:afanoa/SIMPLE_hierch_flex.git SIMPLE && cd SIMPLE && git checkout hetero-tree-wip`.
   The fork has two working branches: `hetero-tree`, the gated tree that is known to work, and `hetero-tree-wip`, which adds the
   ports in progress, this guide and `handoff/` (memories, side-session notes, batch2 scripts). `handoff/README.md` indexes it.
2. Restore Claude's memory (section 10), then read `handoff/memory/MEMORY.md` and its hubs `hub-hetero-tree-current.md` and
   `hub-hetero-tree-verdicts.md`.
3. On batch2, everything below is still in place: builds (section 2), gate fixtures (section 3), harnesses (section 5), the open
   threads' folders (section 7).
4. Before trying an idea, check section 6 (tested and rejected) and the standing rules in section 9.

The fork is public. Nothing personal was pushed: the resume, the GRFP application memories and the GRFP science notes stayed off
GitHub (the GRFP science notes are in the batch2 tarball, section 1.2). The two unpublished manuscript drafts (7.6) are also only in
that tarball.

## 1. Where things are

### 1.1 Git

Fork `git@github.com:afanoa/SIMPLE_hierch_flex.git` (remote `fork` on the Mac). `origin` = `hael/SIMPLE` (upstream): never push there
without the user's explicit permission.

| Fork branch | Content |
| --- | --- |
| `hetero-tree` (`e1d5fdd4c`) | the tree on master's posterior machinery (merge `ee26d9bbf`, gated 2026-11-06) + verdict notes: the version known to work |
| `hetero-tree-wip` | `hetero-tree` + `d142ea8af` depth-tied tree_wiener (`depth_wiener=yes`, default on; gated) + `043b6bcfc`, `c4aa7b562`, `378254f33`, `89604e8ac` per-site split (`site_split=no` default; gating, 7.1) + this guide + `handoff/` |

Everything else from 2026-11-06/07 lives off GitHub, in one git bundle on batch2:
`/mnt/beegfs/elmlund/afan/mac_archive_20261107/fork_archive_20261107.bundle` (114 MB, incremental on `origin/master` and
`hetero-tree`). It holds, as `refs/remotes/fork/<name>`:
- `wip/state-eval` (`eb58a5455`, `refine=state_eval` + `SIMPLE_STATE_PRIOR`, the harness of 5.2);
- `wip/prob-state-marg` (`8d427a32a`, parked, no gain), `wip/tree-ports-charts` (`480bc9f2f`, parked: seams),
  `wip/prior-parentpose` (`df5a3d266`, the pose-anchor scripts);
- `mac/<branch>`: every Mac branch that had commits not on the fork (36);
- `mac/dirty/<worktree>`: `git stash create` snapshots of uncommitted tracked edits in 30 Mac worktrees; `mac/dirty/SIMPLE-main` is
  the main checkout's uncommitted flex2D work (20 files) on top of `mac/flex2D` (`d1a20eff6`);
- `handoff-2026-11-07` (the earlier handoff branch, including 154 MB of the Mac's untracked notes and code) and hetero-tree's
  earlier tip with the guide's history (`1d474d6e0`).

`handoff/ARCHIVE.md` lists them. List with `git bundle list-heads <bundle>`; restore one into a clone that has `origin/master` and
`hetero-tree` with `git fetch <bundle> refs/remotes/fork/wip/state-eval:refs/heads/state-eval`; apply a worktree snapshot onto its
branch with `git checkout <branch> && git stash apply <snapshot sha>`.

### 1.2 batch2

Log in with `ssh batch2.ncifcrf.gov` (the short alias fails host-key verification). Shorthands used below:
`A=/mnt/beegfs/elmlund/afan`, `H=$A/hetero_tree_20261028`, `M=$H/math_reworks`, `L=$M/latent`,
`PY=/mnt/nasapps/production/python/3.12/bin/python3`.

Builds (`<dir>/install/bin/simple_exec`; all passed the fast gate):

| Build | Source | Use |
| --- | --- | --- |
| `$H/merge_master2` | `ee26d9bbf` (= fork `hetero-tree` code) | production tree for every gate and replay |
| `$H/state_eval` | `wip/state-eval` tip (in the bundle) | `refine=state_eval` runs, the oracle |
| `$H/merge_sm` | state-eval + `prob_state_marg` + `write_states` | parked experiment; anchor replays use its `write_states` |
| `$M/build_ports`, `build_ports1` | ports pieces 1+2, piece 1 | depth_wiener gates (`gates/ports_*`) |
| `$M/build_ports3`, `build_ports4` | pieces 1+3; 1+3 with the site_split fix | site_split gates (`gates/ports3_*`, `ports4_*`, `ports5_*`) |

Data and runs:

- Gate harness: `$M/launch_gates.sh`, `$M/sst_gate.sbatch`, `$M/readout.sbatch`; scorers in `$H/fixtures/scripts/` (`qtoy_leaves.py`,
  `closed_nodes2.py`, `make_view.py`, `arm_readout.py`); runs `$M/gates/<arm>_<dataset>[suffix]`; readouts `$M/gates/<arm>_readout.txt`.
- Fixtures: section 3. Toy generators and data `$A/routes/toy55/` (`sim<name>/src.simple`, `truth_qtoy.txt`).
- 692k RNase P trees (all 692,392 particles): `$A/rnasep_klem/stree_port/full692_{m33,base80,base92}_s1/1_hetero_tree`; the closed pair
  is `node_r1111221`'s children `_s1` / `_s2` in m33.
- Python latent framework: `$L` (`latent_states_v4.py`, `flexprior_estep.py`, `calib_flexprior.py`, `flex_refit2_test.py`, `f2r_*`,
  `seval_*`, `orc_*`, `qflex_leaves.py`); per-node dumps `$L/dumps/<key>/` (`b_ptcl3D.txt` = stage B, `r0_ptcl3D.txt` = recB0),
  plane caches `$L/cache/`.
- Thread folders: `$L/oracle`, `$L/flexsub`, `$L/shared`, `$L/global`, `$L/ab_seed`, `$L/seval_prior_runs`, `$M/shrinkB`, `$M/anchor`,
  `$M/mest`, `$M/ppinfo`.
- Session handoffs: `$H/handoff/<session>/` (VERDICT.md or STATUS.md + scripts); every small script of the above in
  `$H/handoff/handoff_scripts_20261107.tgz` (also unpacked in the fork under `handoff/batch2_scripts/`).
- Mac archive: `$A/mac_archive_20261107/mac_untracked_20261107.tar.gz` (1.93 GB): every untracked file of every Mac worktree (3.6 GB
  in the main checkout: result tables, embeddings, the manuscript drafts and their data), the memory folder and the side-session
  notes. File list:
  `handoff/mac_untracked_list.txt`. Extract one path with `tar xzf <tgz> --exclude='._*' "SIMPLE/<path>"` (paths are relative to
  the Mac home; macOS tar stored a `._` metadata companion beside each file, and GNU tar warns about Apple xattr headers: both are
  harmless). Verified on batch2 (job 55013691): the archive lists cleanly, 16,623 files plus their companions.

Rules on batch2: builds on the login node are fine; runs, scorers and tests go through `sbatch` only (csbdevel by default); no
`srun`/`ssh` into compute nodes; no slurm mail. Use the full core cap for tests (80-100 cores across jobs), sized to the job.

## 2. Build and fast gate

```bash
B=$H/<new_build>; mkdir -p $B/build
git clone ~/SIMPLE $B/src && cd $B/src && git fetch https://github.com/afanoa/SIMPLE_hierch_flex.git <branch> && git checkout --detach FETCH_HEAD
cat > $B/build.sh <<EOF
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:\$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
cd $B/build && cmake $B/src -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$B/install -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_ARCHOPT=ON -DUSE_OPENMP=ON > $B/cmake.log 2>&1 && make -j 32 install > $B/make.log 2>&1
echo \$? > $B/build.exit
EOF
nohup bash $B/build.sh > $B/build.out 2>&1 &
```

Then the fast gate, `merge_master2/fast_gate.sbatch` with the paths changed: csbdevel, `-c 8`, `--mem=16G`, 30 min,
`cd $B/src && scripts/run_fast_gate.sh $B/build`, exit 0 expected (15/15 suites on the ports builds).

Traps:
- A backgrounded `&&` chain dies with the ssh session: use a script under `nohup`.
- Every new commit rebuilds everything (~10 min, ~650 objects): `make` rewrites the git-hash header that `simple_defs` includes.
  Iterate on uncommitted edits in the build's `src` (`make -j 32 install` recompiles a few objects), commit in batches.
- A new `.f90` needs `cmake` again (the source list is a configure-time glob).
- Never rebuild an install that running jobs use; make a new build folder per variant.
- Job scripts export the gcc 15.2 `PATH`/`LD_LIBRARY_PATH` themselves (an `sbatch` through plain ssh has no module environment and
  `| tee` hides the failure). Check that a new job's log shows real SIMPLE output within a minute.
- `sbatch --export` splits values at commas: pass keys with commas through the environment (as `launch_gates.sh` does).
- Compare arms only between binaries built the same way. `merge_master2` was built on a transfer node and links OpenMP OpenBLAS
  (`libopenblaso`); the ports builds ran on x6342 nodes with `libopenblas`. Same source, yet the r131 partitions agree only 0.62
  (seed 2) / 0.50 (seed 1) after leaf matching. The porter's `$M/build_ports3.sbatch` (`S=<src> V=<build dir>`: cmake Release,
  `BUILD_TESTS=ON`, 32 cores, ~5 min, runs the fast gate) is the reproducible recipe.

## 3. Gate fixtures: the decisive test for any tree change

```bash
cd $M && ARMKEYS="<the arm's keys>" bash launch_gates.sh <build>/install <ARM> [qtoy2 r131 qtoynull olenull micnull ole qtoyflex]
SEEDX=1 RUNSFX=_s1 ARMKEYS="..." bash launch_gates.sh <build>/install <ARM> r131    # the second r131 seed
sbatch --export=ALL,ARM=<ARM> $M/readout.sbatch                                     # -> $M/gates/<ARM>_readout.txt
$PY -I $L/qflex_leaves.py $M/gates/<ARM>_qtoyflex/1_hetero_tree                     # qtoyflex truth readout (inside sbatch)
```

The launcher exports `SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes` (`sst_gate.sbatch` refuses to start without the split
contrast) and adds `BASEKEYS` (default `oneflex=yes sta0_npeaks=2 junk_screen=no`). Each run writes `run_manifest.txt` (keys, build
sha256, seed, job) and `split_env.txt`. Dataset keys (all with `site_layer=yes site_states=yes`):

| Dataset | Input | Keys | maxdepth / nodepar / seed |
| --- | --- | --- | --- |
| qtoy2 | `$A/routes/toy55/simqtoy2/src.simple` (60k) | `nparts=4 mskdiam=232 minpop=900 kway_maxdepth=2` | 1 / 4 / 1 |
| qtoyflex | `$A/routes/toy55/simqtoyflex/src.simple` (60k) | as qtoy2 | 2 / 4 / 1 |
| r131 | `$H/ancestral/fixtures/r131` (44,334 RNase P) | `nparts=12 mskdiam=232 minpop=3000 kway_maxdepth=-1 split_ntests=460` | 1 / 1 / 2 |
| qtoynull | `$A/routes/toy55/simqtoynull/src.simple` (60k) | `nparts=4 mskdiam=232 minpop=900` | 0 / 1 / 1 |
| olenull | `$A/hetero_tree_regression_20261026/null_controls/ole_null_v1/src.simple` (80.6k) | `nparts=5 pgrp=c2 mskdiam=240 minpop=3000` | 0 / 1 / 1 |
| micnull | `$H/fixtures/pfcrt_micelle/micnull/src_in_fp.simple` (56.8k) | `nparts=5 pgrp=c1 mskdiam=200 minpop=3000 kway_maxdepth=2` | 0 / 4 / 1 |
| ole | `$A/ole_12707/5_refine3D_auto/ole_12707.simple` (80.6k) | `nparts=5 pgrp=c2 mskdiam=240 minpop=3000` | 1 / 4 / 1 |

Pass rules and the reference arm `mm2` (merge_master2 at `ee26d9bbf`, 2026-11-06):

| Fixture | Tests | Pass | mm2 |
| --- | --- | --- | --- |
| qtoynull, olenull, micnull | false splits | one leaf (null acceptance = majority leaf) | one leaf each |
| qtoy2 (6 arm conformers, closed 18, holo 101, junk) | discrete states | NMI(leaves, structure combo) inside the 0.71-0.87 spread of earlier arms (judge against the spread, not one run); closed and holo isolated | 8 leaves, NMI 0.780; best-leaf F1 closed 0.89, holo 0.94, arm combos 0.69-0.94 |
| r131 (~7 % closed reference rows) | HARD GATE: the RNase P closed state is kept | a stage-B map at the closed corner (closed arm >= 0.95, open arm <= 0.6), read on both seeds | seed 2: r2 1.03/0.52 SEPARATES; seed 1 (read on `ports_s1`): r1 1.00/0.53 SEPARATES |
| ole | real-data structure | no collapse | 8 leaves |
| qtoyflex (positive control, section 4) | continuum vs discrete, minority isolation | no rule yet: production fails it | mm2 ran at tree maxdepth 1: best-leaf F1 closed 0.39 (purity 0.27), holo 0.62 (purity 0.47), NMI(leaves, kind) 0.255, the arm cut in two along its coordinate; at maxdepth 2 (the launcher's key; `ports4_qtoyflex`) closed 0.47 (purity 0.39), holo 0.62 |

Arms on record that fail: `mb1` (`prob_beta=1`: r131 never reaches the corner, best 0.94/0.46) and `pb15` (`prob_beta=0.15`: the
corner on 1 of 2 seeds, qtoy2's conformer 1/2 pair merged, NMI 0.748). The r131 corner sits at the readout threshold in every arm
(closed arm 0.91-1.03), and the split test at the corner nodes passes on one or two draws in five, so never call r131 on one tree.

Settled (2026-11-07): the depth_wiener builds (`ports`, `ports3`) give r131 seed 2 a different root split from mm2 (22,484/21,850
vs 23,168/21,166; corner in r1 at 0.96/0.54) while qtoy2 and the nulls are identical. `gates/mm2b_r131` (the mm2 binary again) is
byte-identical to mm2: the tree reproduces run to run. `gates/b0_r131` (`ee26d9bbf` without piece 1, built with the ports recipe) is
byte-identical to `ports`: piece 1 changes nothing before the final maps. mm2 and b0 share the source yet agree only 0.62 after leaf
matching, so the r131 root split is bound to the binary (section 2); a last build separates the BLAS variant from the compile host
(`$H/handoff/ports/STATUS.md`). Judge new arms against a reference built by the same recipe.

## 4. qtoyflex, the positive control

qtoy2 with its six arm stops replaced by a 30-map arm continuum (linear blends of the bracketing deposited conformers), plus the
discrete closed (6 %, conformer 18) and holo (conformer 101) states and junk.

- Generator: `qtoyflex_maps.py` (30 maps), `run_qtoyflex.sbatch` (maps -> `make_oris55` -> drop odd class counts ->
  `QT=qtoyflex bash run_qtoy_sim.sbatch`), in `$A/routes/toy55/` (copies in `handoff/batch2_scripts/routes/toy55/`).
- Data `$A/routes/toy55/simqtoyflex/`; truth `truth_qtoy.txt` (columns `row combo junk A_domain B_blob C_closed jtype conformer s1
  closed_arm holo bfac snrq`; `s1` is the arm coordinate); kinds `$L/truth_qtoyflex_kind.txt`.
- Traps: write the CSVs with `lineterminator='\n'` (CRLF put `\r` into bash paths); `simulate_particles` refuses a class with an odd
  particle count.
- Readings so far: the stage-B root states hold closed at F1 0.17 (recall 0.73, purity 0.10) and holo at 0.62; the finished tree
  reaches closed 0.39 at maxdepth 1 and 0.47 at maxdepth 2. An oracle factorial (7.3) says references, not poses, limit closed here.

## 5. Harnesses: testing a change without a full tree

### 5.1 Node replays with the tree's own calls

`$L/f2r_replay*.sbatch` replays a gate node from its saved step calls (`<step>/hetero_step.sh`: seeds reconstruct3D, stage A, recB,
stage B, state_split), swapping only what the arm changes; at production's seed it reproduces production bit for bit. The replay must
export the run's split environment (`SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes`, the gate's `split_env.txt`): without it
state_split's choose rule flips to the EM groups (qtoyflex: ~10k particles move, holo F1 0.62 -> 0.49). Variants:

- `$M/anchor/pa_replay.sbatch`: extra keys on stage A / stage B (LIST `KEY NODE ARM SEED TRUTH XA XB`, keys joined by `+`).
- `$M/shrinkB/shrinkB_run.sbatch`: stage B as three `maxits=1` calls with `startit=1..3` in one folder (reproduces production: 0.999
  of labels), references edited between calls.
- `$L/flexsub/flexsub_replay.sbatch`: a relocated particle stack (`update_project ptcl_old_root= ptcl_new_root=`).
- Readouts: `$L/f2r_ab_readout.py` (labels and split children vs truth and production); r131 maps `make_view.py` + `arm_readout.py`.

### 5.2 refine=state_eval (bundle ref `wip/state-eval`, build `$H/state_eval`)

Stage B's own refine3D call at the stored poses: each particle's state is the best state reference at its pose; only `state`,
`mi_state` and `corr` are written. It reproduces stage B (qtoy2: 99 % of labels, NMI 0.656 -> 0.653, no drift over 4 iterations;
about a minute per run). With `SIMPLE_STATE_PRIOR=<file>` (rows `iptcl lp_1 .. lp_K`) and `SIMPLE_STATE_PRIOR_ALPHA` (default 1) the
state maximises `alpha*lp_k - T*d_k` (Gaussian log-likelihood at the tables' temperature `T = pi*wsqsum/(2*pftsz)`). Scripts:
`$L/seval_check.sbatch` (LIST `KEY NODE MAXITS [TRUTH]`), `seval_prun.sbatch` (`KEY NODE MAXITS ALPHA TRUTHSPEC`), `seval_prior.py`
(a prior from flex + site factors), `seval_score.py`, `seval_compare.py`. The call is derived from the node's
`stageB/hetero_step.sh` (refine -> state_eval; vol*, prob_neigh_mode, prob_beta dropped; stage B's final maps as vol1..K).

### 5.3 Python latent framework (`$L`; copies in `handoff/batch2_scripts/hetero_tree_20261028/math_reworks/latent/`)

`latent_states_v4.py` (`node_select`, `load_truth`, `readout`, `site_factors`, `xd_fit`/`xd_estep`), `flexprior_estep.py`,
`calib_flexprior.py`, `flex_refit2_test.py` (threaded E-step). Run inside sbatch with `$PY`. Traps: a Python PPCA refit is not a flex
proxy (it fits noise; prior scale ~1e31 on r131): rerun the tree's own fflex call (`$L/flexB_run.sbatch`, `flexB2_run.sbatch`);
`prep_crop` reuses any plane cache whose shape matches.

### 5.4 Real flex_pca on a node at chosen poses

flex_pca refuses a project with several state labels: make one population with
`simple_exec prg=selection oritype=ptcl3D infile=<one 0/1 per row> mkdir=no`. `selection` refuses node-local projects whose stack
segment is the parent's (node_proj=yes); `write_states` (bundle ref `wip/prob-state-marg`, build `merge_sm`) writes states tree-style instead.

### 5.5 Truth readouts

qtoy2 `qtoy_leaves.py <arm>=<tree dir>`; qtoyflex `qflex_leaves.py <tree dir>`; r131 `closed_nodes2.py <outdir> <arm>=<tree dir>`
(8,720 closed reference rows; corner >= 0.95 / <= 0.6). Compare top-N map readouts only at matched N.

## 6. What has been tested (do not retry these forms)

| Idea | Test | Result | Verdict |
| --- | --- | --- | --- |
| flex prior on the refinement's state choice | state_eval + `SIMPLE_STATE_PRIOR`, alpha 1/2 | qtoy2 NMI -0.001..-0.009; r131 closed enrichment -0.01..-0.02; toy holo F1 +0.02-0.03; nulls lose 3-12 % of labels | no |
| same, calibrated by held-out likelihood | Python, frozen poses | gain only with site factors; halves on the full node; switches itself off on r131 | no |
| latent seeding of stage A | node replays | seeds better (qtoy2 0.571 vs 0.309), finished stage B the same (0.650 vs 0.656) or worse (r131 1.02 vs 1.25) | no |
| state label by evidence summed over poses (`prob_state_marg`) | merge_sm build | r131 enrichment 1.450 -> 1.466 (noise), nulls unchanged | no |
| refit flex's class model from the refinement's posterior | Python | olenull flex weight 0 -> 1 (self-reinforcing) | no |
| reliability weights in flex's fit | Python | ~1 everywhere, anti-correlated with particle SNR | no |
| flex at the refined poses (real flex_pca) | flexB runs | own-latent ceiling 0.477 -> 0.504, nothing downstream | no |
| M-estimation, soft weights, tree priors, Student-t, predictive E-step | mest session | nothing beats hard labels + tree_wiener | no |
| composed global chart coordinate | f2r_global | seams between children's charts; no minority isolation; would merge the 692k closed pair | no |
| seam-free tree-wide coordinate (one basis, every particle embedded) | f2r_shared, gate trees | seams and false valleys gone; r131 closed AUROC 0.69-0.72 (d' <= 0.8); not reproducible at final poses; merge rule merges every r131 pair | failing so far (692k pending) |
| posterior temperature | gates at 0.15 / 0.3 / 1.0 | 0.15 and 1.0 lose the r131 corner | keep 0.3 |
| stage-B reference shrinkage toward the node mean (shrinkB: emp2 cross-half Wiener on the state deviations between three single-iteration stage-B calls) | continuation replays, four doses | toys unchanged (qtoy2 NMI 0.656 -> 0.648); on both r131 seeds the split merges the closed corner (full dose, and half dose on seed 1): the shrinkage suppresses the 14.6-10.2 A band the split tests | no: fails the hard gate |
| pose anchor (child poses held near the parent's) | key-only replays, 2 seeds | removes the drift tail (p99 62 -> 9 deg); r131 corner kept, enrichment within seed noise; qtoyflex NMI -0.007; qtoy2 children purity -0.01..-0.11, NMI -0.03..-0.14: it locks in the parent's tempered draw noise | no |
| flex-corrected reconstruction (flexsub: each particle's within-state latent deformation subtracted before recB -> stage B -> split) | replays of the qtoyflex, qtoy2 and qtoynull roots | labels move 0.4-1.1 % (a seed change: 2.9-3.9 %); truth F1/NMI, split and state-map FSCs unchanged; null one leaf. The within-state deformation is ~0.05 % of a particle's flex-band power, and a linear latent cannot un-blur a state map | no |
| better poses vs better references (stage-B oracle) | state_eval factorial, true vs stage-B poses x stage-B vs truth maps | poses are not the limit (true poses with stage-B maps: root accuracy 0.951 -> 0.947); references are: one per truth class lifts the qtoy2 root 0.674 -> 0.919 (ceiling 0.936), r1 0.756 -> 0.839; qtoyflex closed F1 0.17 -> 0.455 | headroom = references |
| depth-tied tree_wiener (port) | 692k + gates + truth maps | nominal leaf FSC 4.7 -> 4.3 A at 692k, state readouts unchanged, gates pass; against the truth maps of the qtoy2 and qtoyflex leaves depth2 never beats emp2 (qtoy2 relative error +0.010 at 20-14 A, +0.016 at 14-10 A; qtoyflex equal within 0.005) | no measurable benefit: the gain is nominal |
| the tree-wide basis as a site dictionary (varimax / ICA localized modes) | six gate trees, input poses | modes localize (0.6-0.8 of a mode's energy in one site) and reproduce, but one mode holds a slice of a site's state: qtoy2 closed AUROC 0.69-0.74 vs 0.92 for the 30-mode span; r131 no mode reaches 0.60 (span 0.63-0.65, site factors 0.65-0.69), and the mode inside the closed-arm mask reads the arm swing | no: discovery and ranking cannot become a lookup |

Findings that frame any next idea:

- The headroom is in the references, not in posteriors or poses (stage-B oracle, 7.3): stage B is held back by its state count
  (the qtoy2 root runs K = 5 for 8 classes, so pairs merge), by split fixed points (r1's 80/20 references are a state_eval fixed
  point) and by junk and minority references. The "32 % confidently wrong" at the qtoy2 root is 27.8 % merged classes plus 4.7 %
  true misassignments, and the true misassignments are unconfident, low-SNR and pose-tailed.
- Stage B is nearly seed-invariant: seeds that agree at NMI 0.03-0.38 end in the same partition.
- Per-particle state evidence on real data sits at the nulls' level (r131 gap 1.8 nats vs 1.4-1.6 on the nulls); on the toys the
  stage-B mistakes are confident (qtoy2: 32 % wrong with median posterior 0.000 on the true class). Map-level statistics carry the
  state information; the split tests already use them.
- Information limits: two equal Gaussians are unimodal below 2 sigma; a 7 % minority needs ~3 sigma for a second mode; the density
  gap test is blind below d' ~2.5 (closed and holo on qtoyflex sit at d' 2.1-2.3, the r131 closed split at 2.2). The supervised AUROC
  ceiling for r131 closed in the flex + site-factor latent is ~0.69.
- Site factors (recB0 state deviations inside site masks) are the only refinement -> flex signal: qtoy2 truth ceiling NMI
  0.48 -> 0.62, r131 closed AUROC 0.65 -> 0.69.
- Child poses drift by one parent iteration's jitter (direction median 2.9 deg) and never accumulate with depth; `prob_neigh_mode=geom`
  confines a direction to the coarse cell of its incoming direction.
- The integrated-node prototype (`$A/paper_inode`) reproduces the production node but integration adds nothing over the stage-B
  analog; focusing references on the node mask is the lever.
- M-estimation has two regimes: toy nodes have one-hot posteriors (soft = hard), real nodes have soft but uninformative posteriors
  (soft dilutes the closed arm). Any weighting by (own leaf, target leaf) is a frequency-blind special case of tree_wiener.

## 7. Open threads and how to resume each

Session handoffs land in `$H/handoff/<session>/`; check there first, it may be newer than this guide.

### 7.1 Ports (fork `hetero-tree-wip`; handoff `$H/handoff/ports/`)

The porter's `STATUS.md` there has exact commands for every step below.

- Piece 1, depth_wiener (`d142ea8af`, default on): GATED 2026-11-07 on `ports_*` / `ports_s1`: nulls one leaf, qtoy2 identical
  (NMI 0.780), r131 corner kept on both seeds (0.96/0.54, 1.00/0.53), final maps deliver depth2. The binary check of section 3
  cleared it (nothing before the final maps changes); the user decides whether it ships on.
- Piece 3, site_split (`043b6bcfc`, default off): an add-only per-site two-state EM on a node about to become a leaf (site-masked
  references, cross-half GAPw with a 3-D hull; `simple_site_split.f90`, the `site_split` program, unit suite; design note
  `doc/implementation_notes/hetero_tree_site_split_2026_11.md` on that branch). Gate 0 (site_split off) = `ports3_*`, identical to
  `ports`. Gate 1 (the worker alone on 17 saved gate nodes, `$M/mest/xs_gate1/<key>/summary.txt`, compared with
  `$M/mest/xs_gate1_pytable.txt` by `xs_gate1_compare.py`; expected 0 DISCRETE except igmfix_root site 6, GAPw 0.09). Run 1
  collapsed every site (fixed in `043b6bcfc`); run 2 found translation splits (group maps docked onto the other half's consensus;
  fixed in `c4aa7b562`); run 3 (build_ports5) reads like the Python on the RNase P leaf nodes, but root nodes keep a strong
  reproducible one-site factor (z1 14-38 where the Python read <= 5.8; all stopped by G2/B, no DISCRETE). build_ports6 logs the
  per-E-step score gap to compare with the Python's (`xs_gate1_gap/`). The first qtoyflex tree with `site_split=yes` ran away (a
  returned node re-filed every poll, 147 leaf lines): fixed in `378254f33` (in build_ports6). Run 3 also missed the IgM positive
  control because the call ran on the EM's coarse grid: `89604e8ac` moves it to the 5.11 A grid its constants were set on (build_ports7,
  run 4 chained). Next: rerun qtoyflex with `site_split=yes` and score it with `qflex_leaves.py`.
- Piece 2, chart outputs: parked (bundle ref `wip/tree-ports-charts`); drop unless 7.2 passes.

### 7.2 Seam-free global coordinate (f2r session; full record and resume commands in `$H/handoff/f2r-global/VERDICT.md`)

One tree-wide basis (all charted nodes' flex modes at 15 A, orthonormalised, top q = 30; `--basis root` as the control), every
particle embedded at its leaf's stage-B pose (`f2r_leafposes.py`) or at its input pose. Jobs 55008359 (six gate trees), 55008360 and
55008361 (692k); outputs `$L/shared/<key>_<all|root>_<final|input>/{report.txt,summary.json}`. Port only if it isolates closed better
than the root chart (AUROC, d', the closed pair's rank among sibling pairs), reproduces across seeds and recipes at least as well,
and supports merges that keep the 692k closed pair distinct.

Gate trees (done 2026-11-06 23:55): the seams are gone, and so are the false valleys on real data and on the continuum. On r131 it
stays weak: closed AUROC 0.69-0.72 (d' <= 0.8; flex's root chart 0.62, the seamed composition 0.67-0.69), and closed-enriched
leaves are no more isolated than other pairs; qtoy2 0.988, qtoyflex 0.843. It reproduces at the input poses (canonical
correlations 0.99-1.00) but not at the final poses (0.87-0.89 on r131: each tree's leaf refinements imprint its own partition).
The Bhattacharyya merge rule would merge every r131 pair and half of qtoy2's distinct pairs. So far it fails the port rule on
the closed state; the 692k trees (the closed pair's rank, cross-recipe CCA) finish the verdict.

### 7.3 What limits stage B: poses, references or noise (oracle session; done; `$H/handoff/stageb-oracle/VERDICT.md`)

A factorial at one node with state_eval (`maxits=1`): poses {stage B, true} x references {stage-B maps, coarse truth at stage B's K,
one map per truth class} (`$L/orc_prep.py`, `orc_run2.sbatch` with env `KEY NODE TRUTH ALLTXT [MUNODE ORC_MINN ORC_MINFRAC]`,
`orc_score.py`, `orc_wrong.py`; runs `$L/oracle/{qflex,qtoy2,qtoy2_r1,qtoy2_r4}`).

- Poses are not the limit. qtoy2's stage-B poses are 2.2-2.6 deg off (median); true poses with stage-B maps move 11-18 % of labels
  but fix as many as they break (root group accuracy 0.951 -> 0.947, r1 0.756 -> 0.758). Poses add 1.5-4 points only once the
  references are pure. Exception: qtoyflex's heavy pose tail, where poses alone give holo F1 +0.065.
- References are the limit. The qtoy2 root runs K = 5 for 8 classes, so pairs merge: one reference per class at the same poses lifts
  class accuracy 0.674 -> 0.919 (ceiling 0.936). r1 (conformers 1 vs 2) sits at an 80/20 reference fixed point: pure references
  0.756 -> 0.839 (ceiling 0.891). r4 is at its ceiling; junk holds closed F1 at 0.896 (0.939 with a junk reference). qtoyflex has no
  closed state (F1 0.17): pure arm/closed/holo references at K = 3 give 0.455, with true poses 0.570.
- So the next work belongs at the split and the reference set (state count, split fixed points, junk and minority references), not
  in posteriors or pose refinement.

Tricks from this session: write true poses or labels into a `.simple` by byte-patching the ptcl3D slots (`orc_patch.py`; header
20 x 5 int64, segment 6, 64 float32 per record; e1 8, e2 9, e3 10, state 23, x 27, y 30); rebuild maps through refine3D's own
PCG/focus/ML pipeline from fixed labels with state_eval and `SIMPLE_STATE_PRIOR` set to 0 on the label and -1e4 elsewhere,
`maxits=2`; child-node frames rotate 1.1-1.3 deg from the truth, so remove the global offset before using true poses.

### 7.4 Refinement priors spawned 2026-11-06 night

- shrinkB, reference shrinkage toward relatives (`$M/shrinkB`, handoff `prior-treeshrink`): closed, worse (section 6). The half
  dose and one-part arms on the toys were still running at handoff (chains 55011976 / 55011977; `sbatch shrinkB_readout.sbatch` puts
  every arm side by side).
- flexsub, removing within-state deformation before reconstruction (`$L/flexsub`, handoff `prior-flexcorr`): at the node's stage-B
  poses, the latent model (flex mean + 10 modes cut to 15 A, per-particle amplitude) gives each particle's deviation from its class
  mean; it is projected (2x padded trilinear slices, gridding-corrected, CTF and amplitude applied) and subtracted from the particle
  image, then recB, stage B and split are replayed on the corrected stack. Arms `full` and `perp` (the subtracted deviation has the
  between-state discriminant directions projected out). Exactness checked (re-read planes = X - D to 5e-8). Closed: no gain
  (section 6; `$H/handoff/prior-flexcorr/VERDICT.md`).
- Pose anchor (`$M/anchor`, handoff `prior-parentpose`): closed (section 6).

### 7.5 Open idea (user, 2026-11-07): one joint coordinate over a node's site flexes, for the data-derived split

Question: "if we have multiple flexes couldn't we have a global coordinate for those (at a particular node)?" Yes, and it avoids the
failure of the composed coordinate. That one failed because tree nodes hold disjoint particles, so each child's chart reaches the
root through its own affine map and the joins are seams. At one node every site chart is fit on the same particles at the same
poses, so a particle has a coordinate in every chart and stacking them `z_i = (z_i^site1, ..., z_i^siteS, z_i^fflex)` needs no
alignment.

What exists already: with `site_layer=yes` (gate keys) each node discovers sites (SVD site factors in site masks) and labels its
particles per site, add-only beside the tree. The layer's 2026-10-21 prototype (`tree_sites.py`) fitted a focused flex per site (the
site mask as the fit mask) and ran a two-class site EM from both tails of the site's leading coordinate: those are the "multiple
flexes" of one node. Earlier results point the same way: on the multi-motion fixture `mm3q` (three independent rigid domains,
8 combinations) the tree fragments the product (1 of 8 combinations is a clean leaf, F1 0.65) while the three root site poles
multiplied give combination F1 0.81 (8 of 8 >= 0.7); on IgM the per-site flex recovers each arm at oracle level, but multiplying many
site labels collapses into the leaves, so the output must stay factorial. Focus on a mask is what lifts minority purity (integrated
node prototype).

What a joint node coordinate could decide that the tree cannot:
1. Independent vs coupled sites: cross-half reproducible dependence (canonical correlations or mutual information) between site
   coordinates. Independent sites are factors, reported per site and never multiplied into leaves; coupled sites form joint states
   worth a split.
2. The split proposal itself: a mixture in the joint coordinate (K by cross-half held-out likelihood) in place of site_em's seeds, with
   the map-level cross-half test still deciding (per-particle evidence alone is at null level on real data).
3. Discrete vs continuous per site: a density gap along each site axis (the charts session's 1-D test with a uniform-continuum null),
   usable where the site chart's d' >= ~2.5. Focused charts should raise d' for a local minority over the whole-node chart; that is the
   bet to test.

Test plan (Python first, existing outputs, no new Fortran):
- Nodes: mm3q root (three independent motions), qtoy2 root (arm + closed + holo), qtoyflex root and r1 (continuum arm + discrete
  minorities), the IgM fixture root, and qtoynull and olenull roots. Check what the gate trees kept: finished trees prune site masks and
  inputs (`SIMPLE_TREE_KEEPMAPS=yes` keeps maps); where pruned, rerun the node's site layer from its `n.simple` (the 2026-10-21
  `tree_sites.py` lineage: `$A/routes/toy55/` and the Mac tarball).
- Per node: each site's flex coordinates with their per-particle precision; dependence matrix between sites; joint mixture; per-site
  gap.
- Readouts against truth: mm3q combination F1 (target >= 0.81, the pole product) and whether independence is recognised; qtoy2 and
  qtoyflex closed and holo F1 (targets: the gate leaves, 0.89/0.94 and, at maxdepth 2, 0.47/0.62); nulls: K = 1 and no dependence.
- Port only if the nulls stay at one state, mm3q's factors are found without fragmenting, and minorities are isolated at least as well
  as the current leaves. On RNase P use only data-discovered sites (as the site layer does): the user has ruled out predefined region
  or multibody bases for this molecule, and the r131 closed corner remains a hard gate.

### 7.6 The paper

3D tree methods paper (3DVA and RECOVAR as the references; RNase P the headline): draft `doc/manuscripts/hetero_tree_2026_11/main.tex`
with `refs.bib` and figures, untracked in the Mac worktree `magical-babbage-701c7b` and kept off the public fork. Get it from the
batch2 tarball: `tar xzf $A/mac_archive_20261107/mac_untracked_20261107.tar.gz SIMPLE/.claude/worktrees/magical-babbage-701c7b/doc/manuscripts`.
Drafted: abstract, introduction, results 2.1-2.4 and 2.7, methods. Pending: 2.5 (integrated node), 2.6 (3DVA/RECOVAR head-to-head),
discussion. The integration claim stands on nothing yet (section 6): the node that works is the production node. Compile with
tectonic (`tectonic -X compile main.tex`). The 2D (flex2D) paper has its own draft, `doc/manuscripts/flex2D_cls_split_2026_10_03/`
in the main checkout, also only in the tarball (`tar xzf ... SIMPLE/doc/manuscripts`).

### 7.7 Next (user, 2026-11-07): rework the site, mask and ranking layer on a test the tree fails

Direction: site ranking and masks are the weakest part of the pipeline. Prove the current tree fails on a harder test, then
explore the flex prior, depth wiener and the global coordinate there, judged by tangible benefit. The present formulation throws
information away at every hand-off: flex latent -> hard cheap children -> their variance -> a 0.2 x max mask -> sites -> two tails
-> poles -> seeds. The oracle (7.3) places stage B's losses exactly in those seeds and references.

The hard test exists: the IgM fixture, `/mnt/beegfs/elmlund/afan/igm_cryolo/tree/fixture_igm/` (`fix/src.simple`,
`null/src.simple`, `truth.txt`, `arm_masks.mrc`, `scripts/score_fx.py`; report
`handoff/side_session_notes/local_5bd42cba-9424-4610-b924-d1ce032a5093/tree_igm_fixture_report.md`). It is the real IgM core
(100k particles, box 200, 3.27 A, real noise) with five 50 kDa Fab-like arms on independent discrete hinges (72 combinations); the
null twin moves the same arms continuously. What the current tree does with it:
- one leaf: sites only in the core, the root mask is flex's solvent imprint, and the site layer finds no arm pole even on injected
  arm sites;
- an arm-focused root finds a1 (NMI 0.13), but the sphere-wide effect floor merges it;
- the best configuration so far (islands + cross-validated site stage B + GAPw, the basis of `site_split`) finds a1 alone (NMI
  0.146), with the null held.

A 10 % arm variant at real Fab proportions (`/mnt/beegfs/elmlund/afan/igm_cryolo/tree/fixture_igm10`, built 2026-11-07, noise
recalibrated) separates an SNR limit from a method limit; its stock roots and status: `$H/handoff/igm10/STATUS.md`.

Second tier: `mm3q` (three independent rigid domains, 8 combinations: the tree has 1 of 8 clean, the site-pole product 8 of 8) and
qtoyflex (closed F1 0.47 at maxdepth 2). The gate fixtures of section 3 stay the guards.

Where the latent work fits:
1. Several flexes per node, one joint coordinate (7.5): a focused flex per discovered site (`pcg_mskfile` = the site mask) at the
   node's poses gives every particle a site coordinate with its precision, stacked across sites with no seams. The whole-molecule
   root latent sees nothing of a 50 kDa arm (arm R^2 <= 0.024, even with an envelope mask); a focused fit is untested on this
   fixture (on the October IgM-like fixture per-site flex recovered every arm).
2. Ranking: one calibrated, comparable number per site (held-out nats per particle of a site mixture against one Gaussian, or GAPw
   on the site coordinate) in place of site_rank's ~45 EM fits and pole margins, whose order alone decides the r131 corner.
3. Flex prior: the site coordinate's posterior as the prior on the site state in a site-focused stage B (harness 5.2). It failed
   where stage B already had the information; weak per-particle evidence, as on the IgM arms, is the regime it was meant for.
4. Global coordinate and depth wiener: the tree-wide basis (7.2) as a site dictionary shared by every node; depth-tied pooling for
   the per-site state maps. Both tested 2026-11-07 and closed (section 6): a site's state is spread over several directions inside
   the site, which focused fits and site factors capture and single localized global modes do not; depth2 has no truth-map gain.

Tests, in order:
- T1 (decisive, ~1 h of cluster): focused flex_pca with `pcg_mskfile` = each arm island of
  `proto/fix_root_stock_recB0_bgyes/sites` (and `arm_masks.mrc` as the oracle-mask control) on `fix` and `null` at the fixture root's
  poses, K = 1, with the tree's own fflex keys. Pass: a site coordinate tracks its arm's truth state beyond the cross-validated site
  stage B (a1 NMI 0.146) and finds a second arm; on the null twin it shows no discrete structure (GAPw no call).
- T2: the joint coordinate's factor-or-coupled calls and per-site discreteness on IgM, mm3q and qtoyflex.
- T3: the new site ranking at the r131 root keeps the corner's site first (the hard gate), qtoy2 unchanged.
- Port only after T1-T3, then the full gate battery.

## 8. Decisions waiting for the user

1. depth_wiener on by default: gated and cleared, but no truth-map benefit (section 6), so the measured case is for `depth_wiener=no`.
2. Drop the chart-coordinate port (bundle ref `wip/tree-ports-charts`) unless the seam-free coordinate passes 7.2's rule (it does not on
   r131 so far).
3. Keep state_eval as a developer harness: add it to `hetero-tree-wip` (after a build and fast gate) or leave it in the bundle.
4. site_split: after gate 1 and qtoyflex.
5. Launch T1 of 7.7 (focused per-site flex on the IgM fixture), which also opens the joint site-coordinate idea of 7.5.

## 9. Standing rules and traps

Science:
- HARD GATE: flex works only at frozen poses; no latent recalculation or Woodbury algebra inside refinement unless rigorously
  proven useful.
- HARD GATE: the RNase P closed state must be retained (the r131 corner, the 692k closed leaves).
- No predefined region or multibody basis for RNase P; the tree is input-free (no seeding from prior knowledge); splits are data
  derived and the landscape must reproduce; science before speed; numerically close is enough, never chase bitwise identity.
- A new coupling must beat section 6's numbers on qtoy2, r131, qtoynull and olenull first.

Cluster and data:
- Never delete cluster data (read, or `mv` aside); never signal SIMPLE processes by name; `scancel` by job id only; never rebuild or
  overwrite scripts under running jobs.
- `simple_exec` keys must follow `prg=`; a failed `chdir` does not abort (the run lands in the project root); finished trees prune
  split inputs.
- Memory sizing from `sacct` MaxRSS; wall time sized to the job.

Git:
- Push only to the fork; commits carry no Co-Authored-By trailer; build and fast-gate on batch2 before pushing code.
- The fork holds two branches: `hetero-tree` (known to work) and `hetero-tree-wip`. Experiments, archives and per-session branches
  stay off GitHub (user, 2026-11-07).

More traps, by topic, in the memory hubs (`handoff/memory/hub-*.md`).

## 10. Restoring Claude's memory and context on a new machine

- Claude Code keeps per-project memory in `~/.claude/projects/<checkout path with every "/" and "." replaced by "-">/memory/`, e.g.
  `/home/me/SIMPLE` -> `~/.claude/projects/-home-me-SIMPLE/memory/`. Copy it in:
  `mkdir -p ~/.claude/projects/-home-me-SIMPLE/memory && cp handoff/memory/*.md ~/.claude/projects/-home-me-SIMPLE/memory/`.
  `MEMORY.md` is the index every session loads; the hubs index the rest.
- Memories name Mac paths (`/Users/atifao/SIMPLE`, `.claude/worktrees/...`) that no longer exist; batch2 paths stay valid.
- `CLAUDE.md` and the repo skills (`.github/skills/`, linked from `.claude/skills/`) are in the repository.
- Side-session status notes: `handoff/side_session_notes/<session id>/`.
- The Mac's untracked files, including the manuscript drafts and result data, are in the batch2 tarball (section 1.2; index
  `handoff/mac_untracked_list.txt`); a 154 MB selection of its notes and code is also on the bundle's `handoff-2026-11-07` ref.
- Not archived: the 1.9 GB of Claude Code session transcripts.
