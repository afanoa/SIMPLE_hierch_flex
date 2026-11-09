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
| `hetero-tree-wip` | `hetero-tree` + `d142ea8af` depth-tied tree_wiener (`depth_wiener=yes`, default on; gated) + `043b6bcfc`, `c4aa7b562`, `378254f33`, `89604e8ac`, `6eb27fa89` per-site split (`site_split=no` default; gating, 7.1) + this guide + `handoff/` |

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
- Compare arms only between binaries that link the same OpenBLAS. The default cmake links OpenMP OpenBLAS (`libopenblaso`, as
  `merge_master2` does); `-DBLA_VENDOR=OpenBLAS` links `libopenblas`. Same source, yet the r131 partitions agree only 0.62
  (seed 2) / 0.50 (seed 1) after leaf matching: the root split there is a near-tie that flips with the library. Compile host and
  code do not matter, and every binary reproduces itself. The porter's `$M/build_ports3.sbatch` (`S=<src> V=<build dir>`: cmake Release,
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
Caution (extremes VERDICT, 7.8): the gate reads each leaf's map from the same particles that selected it, and that imprint is
large. The mm2 corner leaf r2_s2 reads 1.03 / 0.54 on its own data but 0.93 / 0.59 on a fission re-read (random 5 % tails 0.87-0.90,
node mean 0.89-0.90 / 0.59), so no leaf passes the gate honestly. Whether to redefine the gate on fission readouts is decision 9.

Settled (2026-11-07): the depth_wiener builds (`ports`, `ports3`) give r131 seed 2 a different root split from mm2 (22,484/21,850
vs 23,168/21,166; corner in r1 at 0.96/0.54) while qtoy2 and the nulls are identical. `gates/mm2b_r131` (the mm2 binary again) is
byte-identical to mm2: the tree reproduces run to run. `gates/b0_r131` (`ee26d9bbf` without piece 1, built with the ports recipe) is
byte-identical to `ports`: piece 1 changes nothing before the final maps. mm2 and b0 share the source yet agree only 0.62 after leaf
matching. The cause is the OpenBLAS variant (section 2): the base commit built on the same node without `BLA_VENDOR`
(`gates/b1_r131`) reproduces mm2 byte for byte. The closed corner is kept on both seeds under either library. Judge new arms
against a reference that links the same OpenBLAS.

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
| seam-free tree-wide coordinate (one basis, every particle embedded) | f2r_shared, gate trees | seams and false valleys gone; r131 closed AUROC 0.69-0.72 (d' <= 0.8); not reproducible at final poses; merge rule merges every r131 pair; 692k (m33, input poses): the closed pair at d' 0.25, rank 46 of 112 sibling pairs, Bhattacharyya > 0.5 for 109 of 112 | no (the CTF-flip trap, section 9, does not move it: the r131 closed AUROC changes by <= .003 with |CTF|) |
| posterior temperature | gates at 0.15 / 0.3 / 1.0 | 0.15 and 1.0 lose the r131 corner | keep 0.3 |
| stage-B reference shrinkage toward the node mean (shrinkB: emp2 cross-half Wiener on the state deviations between three single-iteration stage-B calls) | continuation replays, four doses | toys unchanged (qtoy2 NMI 0.656 -> 0.648); on both r131 seeds the split merges the closed corner (full dose, and half dose on seed 1): the shrinkage suppresses the 14.6-10.2 A band the split tests | no: fails the hard gate |
| pose anchor (child poses held near the parent's) | key-only replays, 2 seeds | removes the drift tail (p99 62 -> 9 deg); r131 corner kept, enrichment within seed noise; qtoyflex NMI -0.007; qtoy2 children purity -0.01..-0.11, NMI -0.03..-0.14: it locks in the parent's tempered draw noise | no |
| flex-corrected reconstruction (flexsub: each particle's within-state latent deformation subtracted before recB -> stage B -> split) | replays of the qtoyflex, qtoy2 and qtoynull roots | labels move 0.4-1.1 % (a seed change: 2.9-3.9 %); truth F1/NMI, split and state-map FSCs unchanged; null one leaf. The within-state deformation is ~0.05 % of a particle's flex-band power, and a linear latent cannot un-blur a state map | no |
| flex prior on the IgM arms (the focused latent's posterior as the state prior; arm-focused state_eval at the root's poses) | a1, a2, null twin | at most +0.006 NMI at 2-8 points of accuracy (a1 0.168 -> 0.165-0.171; a2 0.103 -> 0.109); it drags oracle references down (a2 0.144 -> 0.116) and on the null overrides the arm objective for 98 % of particles | no: the information is in the references |
| better poses vs better references (stage-B oracle) | state_eval factorial, true vs stage-B poses x stage-B vs truth maps | poses are not the limit (true poses with stage-B maps: root accuracy 0.951 -> 0.947); references are: one per truth class lifts the qtoy2 root 0.674 -> 0.919 (ceiling 0.936), r1 0.756 -> 0.839; qtoyflex closed F1 0.17 -> 0.455 | headroom = references |
| depth-tied tree_wiener (port) | 692k + gates + truth maps | nominal leaf FSC 4.7 -> 4.3 A at 692k, state readouts unchanged, gates pass; against the truth maps of the qtoy2 and qtoyflex leaves depth2 never beats emp2 (qtoy2 relative error +0.010 at 20-14 A, +0.016 at 14-10 A; qtoyflex equal within 0.005) | no measurable benefit: the gain is nominal |
| the tree-wide basis as a site dictionary (varimax / ICA localized modes) | six gate trees, input poses | modes localize (0.6-0.8 of a mode's energy in one site) and reproduce, but one mode holds a slice of a site's state: qtoy2 closed AUROC 0.69-0.74 vs 0.92 for the 30-mode span; r131 no mode reaches 0.60 (span 0.63-0.65, site factors 0.65-0.69), and the mode inside the closed-arm mask reads the arm swing | no (unchanged with |CTF|: best single mode on r131 within .004) |
| stopping rule from the node's flex spectrum (BBP noise edge) and a global r + 1 state bound | 433 nodes of 44 trees + truth (7.8) | nulls carry above-edge pose and contrast directions; the RNase P closed spike sits far below any edge (0.014 at r131); a null-calibrated edge loses qtoy2's close pairs, qtoyflex and r131; r + 1 fails where states exist at r = 0 | no |
| rigid-tangent deflation of the node's flex basis | IgM 5 %/10 %, qtoy2, r131 (7.8) | removes pose but costs real signal (qtoy2 closed .964 -> .944); the IgM false root splits do not follow pose error (AUROC .51-.54). Per-particle profiling keeps the signal (being ported) | no (profiling instead) |
| lower band per site (lever 7: lp 20/25/30 vs 15) | IgM 5 % and 10 %, fix and null (7.8) | fix R^2 falls with lp (10 %: .454 -> .405 at lp 30); the null's coordinates heavy-tail (tied gain .0095 -> .030) | no: keep lp 15 |
| a map-level merge/split test on the data the labels came from (lever 4's merge, swap scoring, view matching) | qtoynull, qtoy2 (7.8 oversplit, hier-ppca) | one-conformer null groups read distinct at z 8-13 (lever 4) and 53-56 (hier-ppca, swap + 300-cell view matching) | no: test on unseen data (band split or fission) |
| the tree as a hierarchical PPCA (B&T node, band-split select, continuum leaves) | qtoy2, qtoyflex, qtoynull, IgM 10 % (7.8) | null 1 leaf; qtoy2 NMI .697 vs tree .780; qtoyflex closed lost (.22 vs .39-.47); continuum vs discrete undecidable from maps; IgM no call | no as a model; keep its test rule |

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
- A cross-half map test of hard-assigned groups is valid only on data the assignment never saw: a band above every assignment and
  pose lowpass, or a data-fission copy (7.8 hier-ppca, fission). The production `state_split` tests inside stage B's 10 A
  assignment band and passes the gate nulls through its locality guard and effect floor; the IgM nulls get through.
- On the IgM fixture the site layer is the bottleneck and discovery is its gate: the latent reads small moving parts only through
  arm-shaped supports (T1), data-only masks reach the true mask where the root's islands find the arm (10 % a1, a3), and
  latent-built references assigned per arm then find the arms (five-arm product .251 vs the tree's .042). Sparse PPCA learns all
  five 10 % supports from the data with clean nulls, and its blocks + one masked round make the best chain measured, data-only
  (product .338 at 10 %, .104 at 5 %; true masks .256 / .098; tree .042 / 0). A discreteness test must still accept each split:
  the seeds cut the null's continuous arms in two.

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
  run 4 chained). Run 4 is clean on RNase P and the nulls (0 DISCRETE in 76 sites) but still misses the IgM control (site 6:
  GAPw .33 vs the Python's .09). The prototype's own `dc_gates.py` on the port's maps, and on the port's labels rebuilt with
  reconstruct3D, both reproduce the port's numbers: the call and the maps are faithful and the LABELS fall short. The port's
  site EM settles slowly (40 % of labels still changing at iteration 7, against `xsbg.py`'s 12-16 % by iteration 8-10), leaving a
  weaker contrast (z1 22 vs 38.5, agreement .81); prime suspect is the per-shell global model scale where the prototype fits a
  per-particle amplitude. Running: the E-step amplitude fit (`SIMPLE_SPLIT_CONTRAST=estep`) on the IgM root and two nulls
  (`$M/mest/xs_gate1_v6`, array 55019481). The hard test runs too: full trees with `site_split=yes` on the IgM fixture and its
  null twin (`gates/ports7_igmfix_ss`, `gates/ports7_igmnull_ss`, jobs 55015788/9), scored into `$H/handoff/ports/IGM_RESULTS.md`
  when they end. qtoyflex with the layer on (`gates/ports7_qtoyflex_ss`) finds no discrete site and equals the control's four
  leaves (closed F1 .473, holo .619).
  `6eb27fa89` (fast gate passed in build_ports10; on `hetero-tree-wip`): the site path fits a per-particle amplitude, which
  removes olenull's brightness factors (z 15 -> < 1), and its cross-half EM runs to convergence (0.5 % of assignments changing,
  10-30 iterations). It still misses the IgM control: v8 site 6 GAPw .43, B 7 (v7, 20 fixed iterations, the same) vs the Python's
  .09, B 40; nulls clean (30 sites, 0 DISCRETE). The IgM hard test with `site_split=yes` (build_ports7): the fixture's tree
  splits 50/50 along nothing (arm NMI <= 0.003, a fail), the null twin stays one leaf (pass); scores in
  `$H/handoff/ports/IGM_RESULTS.md`.
- Piece 2, chart outputs: parked (bundle ref `wip/tree-ports-charts`); drop unless 7.2 passes.

### 7.2 Seam-free global coordinate (f2r session; full record and resume commands in `$H/handoff/f2r-global/VERDICT.md`)

**Checked (2026-11-07):** the r131 and 692k numbers here were computed on ctf=flip stacks with the signed CTF (section 9 trap).
Rerun with |CTF| (`$L/absctf/`, job 55032507; addendum in the VERDICT): the r131 closed AUROC moves by <= .003 (seed 2 span
.631 / .699 -> .632 / .699; seed 1 .649 / .714 -> .652 / .715) and the site dictionary's best single mode by <= .004, because the
basis is band-limited at 15 A, where the affected shells hold little energy. The verdict stands; the 692k cache was not rebuilt.

One tree-wide basis (all charted nodes' flex modes at 15 A, orthonormalised, top q = 30; `--basis root` as the control), every
particle embedded at its leaf's stage-B pose (`f2r_leafposes.py`) or at its input pose. Jobs 55008359 (six gate trees), 55008360 and
55008361 (692k); outputs `$L/shared/<key>_<all|root>_<final|input>/{report.txt,summary.json}`. Port only if it isolates closed better
than the root chart (AUROC, d', the closed pair's rank among sibling pairs), reproduces across seeds and recipes at least as well,
and supports merges that keep the 692k closed pair distinct.

Gate trees (done 2026-11-06 23:55): the seams are gone, and so are the false valleys on real data and on the continuum. On r131 it
stays weak: closed AUROC 0.69-0.72 (d' <= 0.8; flex's root chart 0.62, the seamed composition 0.67-0.69), and closed-enriched
leaves are no more isolated than other pairs; qtoy2 0.988, qtoyflex 0.843. It reproduces at the input poses (canonical
correlations 0.99-1.00) but not at the final poses (0.87-0.89 on r131: each tree's leaf refinements imprint its own partition).
The Bhattacharyya merge rule would merge every r131 pair and half of qtoy2's distinct pairs. It fails the port rule, and the 692k tree
confirms it (m33 at the input poses; the cross-recipe job was cancelled as non-decisive): the closed pair sits at d' 0.25, rank
46 of 112 sibling pairs, the merge rule would merge 109 of 112 pairs, and neither the root split nor the pair axis has a valley.
Closed.

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
recalibrated so the particles keep the real band correlation; status, resume and scoring: `$H/handoff/igm10/STATUS.md`). Interim
root results (merge_master2, 03:20) show the current tree failing a second way: the fixture's root splits in two weakly along
a1 (NMI a1 0.072, other arms <= 0.015, groups 65-66 % pure in a1; the guard keeps the split), and the NULL twin's root splits
50/50 along the continuous arm angles (largest leaf 0.505, a null-gate fail; angle eta^2 a1 0.050, a3 0.032, a5 0.032). At 5 %
(build_igm4) both roots were one leaf, the guard removing a nuisance split. The same-build 5 % roots (jobs 55016840/1) decide
whether this is the 10 % signal or the build. They did, and the production binary fails at 5 % too: on merge_master2 the 5 % fix
root splits 50/50 along nothing (arm NMI <= 0.004) and the 5 % NULL splits 64/36 (pre-guard 3 groups; angle eta^2 0), where
build_igm4 held both at one leaf and the site_split build (the other OpenBLAS) holds the null at one: the IgM null root decision
is a near-tie, like r131's split. The root's whole-molecule latent is still blind at 10 % (arm R^2 0.001-0.006,
LDA at chance). T1 at 10 % (true arm masks a1-a5 plus the stock mask) was running at handoff (jobs 55020271 fix, 55020268 null);
results land in `fixture_igm10/t1_sitefocus/runs/<twin>_<mask>_a<k>/score.json`, the table command is in the igm10 STATUS.md.
T1 at 10 % (true masks): every arm read strongly, R^2 0.40-0.56, LDA NMI 0.28-0.44, and now UNSUPERVISED GMM NMI a1 0.30, a2
0.21, a3 0.29, a5 0.22 (a4 0.003), no leakage; the whole-molecule latent reads ~0; the null's focused coordinates track the
continuous angle (R^2 0.25-0.31) without clusters; per-particle discreteness still does not separate fix from null. At real Fab
scale a focused latent with a good mask recovers four of five arms unsupervised, where the tree finds a1 at 0.072.

Second tier: `mm3q` (three independent rigid domains, 8 combinations: the tree has 1 of 8 clean, the site-pole product 8 of 8) and
qtoyflex (closed F1 0.47 at maxdepth 2). The gate fixtures of section 3 stay the guards.

Where the latent work fits:
1. Several flexes per node, one joint coordinate (7.5): a focused flex per discovered site (`pcg_mskfile` = the site mask) at the
   node's poses gives every particle a site coordinate with its precision, stacked across sites with no seams. The whole-molecule
   root latent sees nothing of a 50 kDa arm (arm R^2 <= 0.024, even with an envelope mask); a focused fit is untested on this
   fixture (on the October IgM-like fixture per-site flex recovered every arm).
2. Ranking: one calibrated, comparable number per site (held-out nats per particle of a site mixture against one Gaussian, or GAPw
   on the site coordinate) in place of site_rank's ~45 EM fits and pole margins, whose order alone decides the r131 corner.
3. Flex prior: the site coordinate's posterior as the prior on the site state in a site-focused stage B (harness 5.2). Tested on
   the IgM arms 2026-11-07, the regime it was meant for: no gain (section 6).
4. Global coordinate and depth wiener: the tree-wide basis (7.2) as a site dictionary shared by every node; depth-tied pooling for
   the per-site state maps. Both tested 2026-11-07 and closed (section 6): a site's state is spread over several directions inside
   the site, which focused fits and site factors capture and single localized global modes do not; depth2 has no truth-map gain.

Results so far (2026-11-07, `$H/handoff/sitefocus-gates/VERDICT.md`, `tables.md`):
- Focused per-site flex does not beat the whole-node latent where that latent already sees the state: r131 closed rows, cross-half
  LDA AUROC 0.53-0.60 per site vs 0.61 for the node (both seeds; site 2 covers the closed arm); qtoy2 best site closed 0.94 / holo
  0.93 vs 0.98 / 0.995. Node + all sites is never worse (r131 0.617), and qtoyflex gains (closed 0.786 -> 0.825 with site 1, 0.845
  with all sites). On the IgM fixture, where the whole-molecule latent reads 0.000 for every arm, a fit focused on arm a1's true
  mask reads a1 (T1 below). So focused coordinates belong beside the node's, not instead of them.
- T1, IgM fixture (`/mnt/beegfs/elmlund/afan/igm_cryolo/tree/fixture_igm/t1_sitefocus/T1_VERDICT.md`; generic scripts
  `scripts/sitefocus_run.sbatch`, `sitefocus_score.py`): with the TRUE arm masks a focused fit sees every arm (R^2 0.19-0.31,
  cross-validated LDA NMI 0.11-0.21: a1 0.211, a2-a5 0.11-0.14, against the cross-validated site stage B's a1 0.146, a2 0.047,
  a5 0.039, a3/a4 none; best unsupervised single-axis split a1 0.144, a3 0.119, a5 0.108, a2 0.104, a4 0.006; no leakage between
  arms), while the whole-molecule latent reads 0.000 on every arm. With the tree's own ISLAND masks it mostly fails (R^2
  0.04-0.11, LDA NMI 0.016-0.063): each island holds only 13-17 % of its mass on its arm, the rest is solvent, and the coordinates
  turn heavy-tailed (kurtosis 5-28) in both twins. Per-particle discreteness is unreadable at this SNR (tied-variance gain ~0 for
  the discrete and the continuous arms alike), so the discrete-vs-continuous call stays map-level (GAPw). Verdict: the mask is
  the bottleneck; site discovery must deliver arm-shaped masks before the latent can help.
- T1b, data-only arm masks (`fixture_igm/t1b_masks/T1B_VERDICT.md`, `scripts/t1b_masks.py`): trimming each island to the
  molecule's support (the tree's 0.1 x max rule) doubles the island fit on every arm, means over the five arms R^2 0.066 -> 0.127
  and LDA NMI 0.032 -> 0.070 (true masks 0.231 / 0.142); a1 reaches R^2 0.244, LDA NMI 0.162, past the cross-validated site
  stage B's 0.146. Purity matters more than coverage: the trimmed masks are as arm-pure as the true ones (0.77) but cover a third
  of the arm, because the arm is smeared in the node mean (73-90 % of each true arm mask lies outside that support); a looser
  threshold restores coverage, loses purity and scores worse. Trimming also removes most of the heavy tails (kurtosis 3-8).
  Refining the island by its focused fit's own mode energy (`self`) adds nothing over trimming (mean R^2 0.120 vs 0.127; a1
  0.211 / 0.138 vs 0.244 / 0.162), and `trim+self` masks cover only 20-30 % of the arm (fits running at handoff, job 55019659,
  scoring themselves; table `t1b_masks/scripts/t1b_table.sbatch`). Trim is the best data-only mask; a mask that follows the arm's
  low-density swept volume is still missing, and it is the next step.
- The ranking number is not calibrated yet. With free variances heavy tails buy it (8 of 15 qtoynull sites positive, up to 0.295
  nats/particle on a near-empty site, kurtosis 48; at r131 its fitted tail is not the closed rows, AUROC 0.51). With tied variances
  it is clean on the null (all 15 sites < 0) but ranks site 3 above the closed-arm site 2 at r131 seed 2 (hard gate). The tree's
  own pole rule never keeps a pole on site 2 and, at seed 1, seeds from site 1, which covers neither arm.
- The tree's site choice misses the states that matter on all three fixtures with truth: qtoy2 picks s1 (closed arm = s3 0.72,
  s5 0.60; holo excess = s2), qtoyflex picks s3 (closed arm = s1 0.72; holo = s2), r131 seed 1 picks site 1 (closed arm = site 2,
  0.77-0.81; site 1 covers 0.05); r131 seed 2 seeds from its fflex states (no site gives a split). Both r131 trees still deliver the
  closed corner at r2. Oracle-mask controls (truth masks cut to the site-mask shape): qtoy2 closed 0.974 (d' 2.72) vs the tree's
  best-covering site 0.881 (1.64) vs node 0.976 (2.85): placement costs 0.04-0.09 AUROC and focus gains nothing; qtoyflex closed
  0.835-0.852 vs node 0.757-0.786: focus helps, the tree's site keeps little (0.764); holo masks fall far below the node on both toys.
  On the oracle closed masks the 2-Gaussian tail is the closed state (tail AUROC 0.96 qtoy2, 0.76 qtoyflex); at the tree's r131
  sites it is not (0.51). At r131 the oracle masks change nothing: closed AUROC 0.604 (mask covering the closed-arm region) and
  0.613 (dilated readout region) on z1..10 vs the node's 0.612, the tree's site 2 0.599; placement costs at most 0.03 and focus
  gains nothing, and even there the 2-Gaussian tail is not the closed state (0.50-0.53). The RNase P closed state is weak in any
  per-particle flex coordinate; the r131 corner lives at map level.
- Flex prior on the IgM arms (`$H/handoff/igm-prior/VERDICT.md`, runs `$M/latent/igm_prior/`): no gain (section 6). But references
  built from the focused latent's split and assigned with the arm mask (state_eval, root poses) reach a1 NMI 0.168 / accuracy
  0.714, past the latent split (0.144) and the cross-validated site stage B (0.146); references from the truth labels add 0.03-0.04
  more. So the route is arm-shaped masks -> focused latent -> references -> site-focused assignment, not a per-particle prior.
  Arms a3-a5 not run (resume command in the verdict).
- Traps: flex_pca refuses `fflex/f.simple` (it carries the delivered state labels): use the node's `n.simple`; one compute node
  (cn121) ran flex workers ~100x slow.

Tests, in order:
- T1 (decisive, ~1 h of cluster): focused flex_pca with `pcg_mskfile` = each arm island of
  `proto/fix_root_stock_recB0_bgyes/sites` (and `arm_masks.mrc` as the oracle-mask control) on `fix` and `null` at the fixture root's
  poses, K = 1, with the tree's own fflex keys. Pass: a site coordinate tracks its arm's truth state beyond the cross-validated site
  stage B (a1 NMI 0.146) and finds a second arm; on the null twin it shows no discrete structure (GAPw no call).
- T2: the joint coordinate's factor-or-coupled calls and per-site discreteness on IgM, mm3q and qtoyflex.
- T3: the new site ranking at the r131 root keeps the corner's site first (the hard gate), qtoy2 unchanged.
- Port only after T1-T3, then the full gate battery.

### 7.8 Levers (2026-11-07; the user excluded signal subtraction)

Each writes to `$H/handoff/<name>/` (VERDICT.md, or STATUS.md with job ids and resume commands). State at 21:36 EST:
- `igm-masks` (levers 1 + 7): data-only arm masks from reconstructions binned along the focused coordinate (the arm's swept
  volume, after cryoDRGN's landscape analysis), and lower bands per site (lp 20/25/30 A); IgM 5 % and 10 %.
  VERDICT (`$H/handoff/igm-masks/VERDICT.md`, tables in `results.md`): lever 1 works at 10 %, where the root's islands find the
  arm, and not at 5 %; lever 7 is a no.
  - 10 %, island arms: the noise-calibrated swept mask (`zh`: bin maps > 100 sigma, anchor-connected, data-only) reaches the true
    mask: a1 R^2 / LDA NMI .547 / .433 (true .556 / .444, trim .483 / .365), a3 .474 / .341 (true .455 / .320). On the null twin it
    reads the angle as well as the true mask (R^2 .277 vs .274) and adds no discrete structure (tied gain .010, true-mask level).
  - 10 %, arms without an island (a2, a4, a5; started from transferred masks): iteration 2 lifts a2 .197 -> .419 (true .425), a4
    only .198 -> .226 (true .399), a5 .207 -> .289 (true .434).
  - 5 %: no variant beats the T1b trim on average; only a1 gains (zh1 .296 vs trim .244, true .313). The bin maps carry the flex
    modes' selection bias inside the start mask (1,400-2,400 sigma vs <= 196 outside), so the briefed 0.2 x peak rule re-makes the
    start; where the start fit's z1 is not the arm (a2, a4) the bins build the wrong mask and iteration 2 hurts. No truth-free
    adoption gate found.
  - Discovery is still the gate: the 10 % fix root's recB0 islands cover a1 and a3 only (sparse-ppca below finds all five).
  - Lever 7: no gain anywhere (mean fix R^2 at lp 15/20/25/30: 10 % .454/.445/.436/.405, 5 % .231/.229/.219/.192; null flat);
    lower bands heavy-tail the null coordinates (10 % tied gain .0095 -> .030). Keep lp 15.
  - Masks (under `$A/igm_cryolo/tree/`): `fixture_igm10/igm_masks/masks/zh1_{fix_a1,fix_a3,null_a1..a5}.mrc` and
    `zh2_fix_a{2,4,5}.mrc`; at 5 % the T1b trims stay best except `fixture_igm/igm_masks/masks/zh1_fix_a1.mrc`. The recipe
    (parts in-job, `SIMPLE_SEED=1`, cn121 excluded) reproduces T1 to 3 decimals; the 3-D GMM NMI is recipe-fragile (.000 vs T1's
    .300 on the same coordinates): judge by R^2, LDA and the best axis.
- `igm-refs` (lever 3, the oracle session): seed each arm's states from the arm-focused latent, build references from those
  labels, re-assign with the arm mask (state_eval), iterate. VERDICT (`$H/handoff/igm-refs/VERDICT.md`, 48 tables in `results/`):
  - 10 %, five-arm product vs the 72 combinations, seed -> round 1 -> round 3: data-only zh1 masks .165 -> .225 -> .251, true
    masks .222 -> .235 -> .256, tree .042. Per arm at round 1, zh1 (true): a1 .290 -> .356 (.351), a3 .315 -> .349 (.346), a4
    .000 -> .085 (.013) and .194 by round 3 (.139: the data-only mask beats the dilated true mask); the 3-state arms a2 / a5 seed
    weaker (.113 / .112 vs true .215 / .219) and reach .168 / .172.
  - 5 %: every zh1 seed <= .001; the rounds recover a1 like the true mask (.054 round 1, .131 round 3; true .056 / .137) and
    little else; product .044 vs .098.
  - Null: 10 % product .073 -> .077 -> .093 (true masks .059 -> .068 -> .091); 2-state null arms <= .007 at round 1, 3-state
    +.013-.016; 5 % null .006 -> .012.
  - Rule: one round, 2-state splits are clean; any split needs a map-level test before the tree accepts it, and that test must
    use data the labels never saw (hier-ppca and fission below). The discrete information sits in the focused seed.
  - Part 3, on sparse PPCA's supports (FINAL; tables `results/spf10_*`, `spn10_*`, `spf5_*`; jobs 55032554-56, 8 GB workers):
    mask = each block's support resampled to the 200^3 grid with T1's soft edge; seed = the block's own coordinates (sparse PPCA's
    GMM recipe, reproducing its NMIs exactly). Five-arm product vs the 72 combinations at round 1: 10 % .338 (seed .333, round 3
    .331), true masks .235 (.256 at round 3), zh1 .225 (.251), tree .042; 5 % .104 (seed .050), true masks .075 (.098), zh1
    .025, tree 0, with a5 missed and the non-arm block included. 10 %: the round lowers the 2-state arms slightly toward the
    refinement's fixed point (a1 .394 -> .378, a3 .365 -> .349, a4 .324 -> .318) and raises the 3-state ones (a2 .282 -> .294,
    a5 .295 -> .318). 5 %: the round is decisive (a3 .012 -> .139, a4 .011 -> .117, a2 .028 -> .084, a1 .178 -> .192), and the
    sparse supports beat the dilated true masks as assignment masks. Null (10 %, 8 blocks): sparse PPCA's seeds already split the
    continuous motion (.15-.16 vs angle bins; product .180) although its own held-out test calls it continuous; the rounds add
    0 to -.006 on 2-state high-purity blocks, +.004 / +.007 on low-purity second blocks, +.012-.034 on 3-state blocks (one a5
    block climbs .002 -> .117 by round 3); 5 % null: 0 sites. Rule to port: sparse-PPCA blocks + one masked round as the
    data-only seeding, and a discreteness test (sparse PPCA's held-out gain, GAPw / cross-half, or the fission test) decides
    which splits the tree accepts.
  - Part 4, the per-site discreteness test (FINAL; VERDICT.md "Part 4", per-block files `results/part4/`; jobs 55033299-302,
    55035034): GAPw (`dc_gates.py`, constants unchanged) on the round-1 K = 2 groups of each sparse block. Of 10 should-accept
    blocks (9 IgM arms + qtoyflex closed) and 10 should-reject (8 null blocks, the 5 % non-arm block, qtoyflex's arm block):
    GAPw 9/10 power, 0/10 false accepts; the 2- vs 1-Gaussian held-out gain (bootstrap-calibrated) 8/10 and 4/10; a spike-vs-
    continuum family test 8/10 and 8/10. IgM: all 9 arms DISCRETE (10 % .03-.24, 5 % .05-.27, 3-state arms split in two
    included), null blocks .48-1.56 or no call, the non-arm block no call; the null's true angle-median split is not discrete
    either (.47), so the separation is the test's, not the seed's. Margin .27 vs .47; the .28 threshold is the only calibrated
    number (the 5 % IgM null, so IgM is partly in-sample). qtoyflex fails both ways even with truth groups: closed vs rest is no
    call (lobes 9 A apart, G2 needs 10; the arm continuum fills the corridor), and the arm continuum split at its true median
    reads DISCRETE (.13, lobes 30 A apart: a large curved swing leaves the chord empty). The per-particle tests detect
    non-Gaussianity (Gaussian null -.0025 +- .0002 nats/particle) or are misspecified (real continua prefer 3 spikes; a d' ~2 arm
    fits a blurred segment). A fission map test accepts carved continua by construction (hier-ppca: adjacent arm thirds z 20.5 /
    5.2), so it was not run. Rule: GAPw accepts a site split on IgM-like sites, no call = one leaf; large curved continua need a
    node-specific null, and a minority state next to a continuum needs a test that does not require an empty gap.
  - Part 4b (FINAL; VERDICT.md "Part 4b", `results/part4b/`): GAPw on every knob-free p3 block (round-1 K = 2 groups) accepts
    the five arms' main blocks and nothing else:
    - 10 % fix: 5/5 DISCRETE (a1 .04, a3 .07, a4 .08, a2 .22, a5 .19);
    - 10 % fix extras b5, b7, b8 (a5's main block is b6): 0/3 accepted;
    - 10 % null: 0/13 (.48-1.65);
    - qtoyflex: the closed-reading block b3 gives no call (lobes 3 A; its round-1 split isolates closed no better than p2's,
      NMI .094).
    So GAPw on every block is safe on IgM-like sites, and the extra blocks need no filtering before the chain.
  - Part 5 (FINAL; VERDICT.md "Part 5"; jobs 55036268, 55036321-23): a knob-free saturation test with data fission.
    - Method: fission in memory (tau 1) on the framework's cached Fourier planes; a y1 log-likelihood score orders the particles
      into four groups (each side of the input split halved at its median); y2 maps give cross-half squared steps q_A, q_mid, q_B
      in the site mask. DISCRETE iff z(q_mid) >= 3 and one side has r = q_side / q_mid with r + 2 SE < min(1, (2w)^2), w = that
      side's fraction. A discrete site saturates (its groups are mixtures of two maps); a uniform continuum has equal steps
      (r = 1) whatever its curvature.
    - Results: qtoyflex closed vs rest (truth) DISCRETE (z 4.6, rest side r .36 +- .10; GAPw no call); the arm continuum at
      its true median not discrete (r 1.37 / 1.24; GAPw DISCRETE .13); IgM 10 % fix a1 DISCRETE from truth (.73 / .72, z 23.7)
      and from round-1 labels (.79 / .73); null a1 not discrete (1.39 / 1.15 truth, 1.50 / 1.24 round 1); a random split no call
      (z_mid -0.6; 23.2 without fission, so the imprint is as large as a real middle step).
    - Costs and limits: fission halves the score's information (10 % a1 d' 2.27 -> 1.55, r .50-.55 -> .72-.79; closest accept
      .99 vs 1.00), so the 5 % arms may fall to no call; a mode at one end of a continuum reads as a state; a wrong split
      direction gives no call, not a rescue (qtoyflex's realistic round-1 closed split: the labels do not separate closed); the
      maps are the framework's ML reconstructions, not SIMPLE's PCG maps.
    - In the tree: one seeded fission draw per particle on the fly, the score at frozen poses on y1, four forced-label
      state_eval reconstructions on y2, three steps.
  - Part 5 addendum (FINAL; jobs 55037212, 55037693-94, 55036709, 55036868):
    - On the round-1 labels of the p3 blocks, the saturation test accepts 3/5 main arm blocks, 0/3 extras, 0/13 null blocks and
      0/3 qtoyflex blocks.
    - It loses both 3-state arms split in two (a2 r 1.10 / 1.37, a5 1.52 / .87): the K = 2 boundary runs through the middle
      state, so the steps are collinear and equal, as in a continuum. On qtoyflex the round-1 groups put the holo / closed jump
      inside a side.
    - The end-peaked control (IgM null a1 with a 15x occupancy fall) does not read as a state (dense-side r 1.04 / 1.08). The
      score's ordering noise blurs the gradient; a sharper score would bring the limit back (true-s1 steps .66).
    - Reading: the saturation test fixes GAPw's qtoyflex failures only with true groups. It is a check for curved continua,
      not a replacement for GAPw.
  - Part 6 (FINAL; VERDICT.md "Part 6", `results/part6/`; jobs 55038086-89): both tests on the same round-1 groups.
    - Power on IgM arm blocks: GAPw 14/14, saturation 6/14 (10 % p2 3/5, 5 % p2 0/4, 10 % p3 3/5).
    - False accepts: 0/25 for each (10 % null p2 0/8, 5 % non-arm 0/1, p3 extras 0/3, null p3 0/13).
    - The saturation test loses power in two ways. Fission halves each particle's information, so at 5 % the four groups'
      composition is a ramp (a1: .78 / .57 / .41 / .21) that reads like a continuum. And a K = 2 split of a 3-state arm leaves
      the middle state straddling the boundary.
    - GAPw's margins: arms up to .27, nearest reject .36 (p3 extra b7), null blocks from .47.
    - qtoyflex: neither test calls closed from round-1 labels, so the site layer will not confirm RNase P's closed corner; the
      tree's own splits must carry it.
    - The end-peaked controls give no call or not discrete (qtoyflex subsets are too small to reach z 3).
    - Decision 8: the site layer uses GAPw.
- `oversplit` (lever 4, the f2r session): over-split (K+2 to 2K from the latent, plus a junk component), assign at fixed poses,
  merge by a cross-half map test; qtoy2 root and r1 first, then qtoyflex, r131 (hard gate), qtoynull.
  VERDICT (`$H/handoff/oversplit/VERDICT.md`): the best stage-B variant measured. Over-split (2K seeds from the node's deconvolved
  flex mixture, agglomerated closest-first) + a data-only junk reference (robust z of -corr/-contrast > 3) + one cross-half
  re-assignment (each half labelled with the other half's maps) takes the qtoy2 root from class accuracy .674 to .803 (NMI
  .662 -> .738, closed F1 .55 -> .95, holo kept whole, K 9; oracle .919); qtoyflex gains a closed-enriched state (F1 .17 -> .45).
  The junk reference never costs class accuracy (junk F1 .14 -> .73 at r1) but grows into a 7.7k low-contrast sink on r131. Limits:
  r1 stays at its 80/20 fixed point (.757); r131's closed corner is neither found nor lost; seeds are the limit (seed NMI .53 root,
  .21 r1, .004 r131). The map-level merge test fails the null by self-selection (five states of qtoynull at z 8-13) and tilt-band
  splits look real (z 7). hier-ppca found why and fission fixes it (below): the test has to run on data the labels never saw.
- `pose-tangent` (lever 5): the node mean's rigid-body tangents projected out of the fitted flex basis, particles re-embedded.
  VERDICT (`$H/handoff/pose-tangent/VERDICT.md`): basis deflation fails. The node latent is pose-contaminated (the IgM roots'
  leading axis is ~30 % shift error; the nulls' only mixture structure is pose). Deflating the basis removes it but costs real
  signal (qtoy2 closed .964 -> .944: domain motions carry rigid components); flex2D's plain-metric form does nothing; profiling
  the six tangents per particle (Schur complement) removes pose cleanly and keeps the signal, but the tree's IgM false root splits
  do not follow pose error (AUROC .51-.54): hygiene, not the failure.
  - Port of the per-particle profile (user: "Port it"; pose-tangent2 session; `$H/handoff/pose-tangent/STATUS.md`): patch
    `$H/handoff/pose-tangent/pose_profile_f3eab71ce.patch` on f3eab71ce, in `embed_latents_with_contrast`
    (`src/main/flex/fit/simple_flex_pca_embed.f90`). It was on by default; since `4afe93b80` it is off, and
    `SIMPLE_FLEX_POSE_PROFILE=yes` turns it on (decision 7). Build
    `$M/build_poseprof`, fast gate 15/15. At both IgM 5 % roots flex's own coordinates become pose-free (shift-error R^2 .003 vs
    .35 off); off reproduces the old coordinates (R^2 .998-1.000). Same-binary on/off gate trees (`$M/gates/pp_{on,off}_*`,
    readouts `pp_{on,off}{,_s1}_readout.txt`) PASS: qtoynull, olenull and micnull one leaf with identical partitions; qtoy2 8 leaves
    in both, structure NMI .811 vs .780; the r131 closed corner kept on both seeds with the profile on (SEPARATES(map)); the r131
    root partition moves between arms (NMI .17, the known near-tie), 4 leaves in both; 98.9 vs 96.1 core-h. COMMITTED on
    `hetero-tree-wip` as `5f80720e0`.
  - Site shift profile (profile-sub-pixel-shift session; `$A/igm_cryolo/tree/fixture_igm/xsb_shiftprof/VERDICT.md`, classifier
    `xsb_shiftprof/scripts/xstageb_sp.py`): two shift columns per particle Schur-complemented out of every group score, with the
    masked in-plane gradient of the groups' MEAN reference as the tangent, kill the IgM core-patch shift split (z 26.7 / 25.8 ->
    4.7 / 1.4 swap / independent; label R^2 on the shift error .32 -> .000) and keep a1 discrete (z 38.1 -> 33.8, NMI .146 ->
    .129). A ~1 px prior fails (core z 12): the site's own shift information (~.5-1 px^-2 at the core patch, .03 at an arm) is no
    better than the prior. Untested: a2-a5, the null arms, a Fortran port.
- `sparse-ppca` (lever 6, sparse-ppca2 session; VERDICT YES: `$H/handoff/sparse-ppca/VERDICT.md`, runs `$M/sparse_ppca/runs/`): a flex-style PPCA on the latent framework (node mean, per-shell scale and noise, per-particle amplitude, 12 fixed
  nuisance templates with per-particle coefficients) with a per-shell empirical-Bayes prior; a group lasso over supervoxels
  selects each block's support and the modes are refitted without the penalty (lambda recalibrated every iteration: scale-free);
  sites come from a greedy pursuit over seed balls on the node mean's envelope, by held-out (odd-half) likelihood gain against
  2 SE and 12 solvent controls (after structured sparse PCA, Jenatton, Obozinski & Bach 2010).
  - Validation: at zero penalty it reproduces flex's focused fits within 5-12 %; synthetic (3 planted sites x 2 modes): one block
    per site, 98-99 % of its energy on its site, the full span's R^2 (.64-.66), supports 1.3x the site.
  - qtoyflex: two blocks on one 3,232-voxel arm support; closed / holo AUROC .896 / .926, against oracle-mask fits of the same
    model .852 / .851, global PPCA .817 / .883 and flex's node latent .757 / .861.
  - IgM 10 % fix (`runs/x10fix_p`): all five arms learned from the data alone, one 2-mode block per arm, 95-97 % of a block's
    energy on its arm, 99 % coverage of the arm's state variance, no leakage. R^2 / GMM NMI / LDA NMI: a1 .523 / .394 / .394, a2
    .393 / .282 / .285, a3 .484 / .365 / .372, a4 .441 / .324 / .337, a5 .395 / .295 / .297: the true-mask fits' level, GMM NMI
    above T1 flex's on all five arms.
  - Nulls (after two fixes: one lambda per group in the rotation step, and a leave-one-block-out check, without which the last
    block of a null survives): the 10 % null twin's arms are found as continuous motion with no discreteness (held-out 2- vs
    1-Gaussian gain < 0); the 5 % null and qtoynull give 0 sites. IgM 5 % fix: 4 of 5 arms at T1-level readouts (a5 missed) plus
    one weak non-arm block that passes the leave-one-out check by a hair (1.3 nats vs controls 0.7) and reads nothing.
  - Why it works (BBP): an arm's whitened eigenvalue (~2.3) is below a global PCA's noise edge (~4.8) but above an arm-sized
    support's (~1.0), so global latents and rotations cannot see arms; local seeds plus sparsity can.
  - Knobs (VERDICT table): kappa 3 (in units of the group's pure-noise gradient), 3^3-voxel groups, seed balls, 2 modes per site,
    iteration counts, and `pursuit_rel` 2 %, a safety stop that ended the IgM 10 % pursuit at round 6 (final run `x10fix_p2`: 196
    nats vs a cut of 203, solvent controls 0.1) before a sixth block was tried.
  - Knob-free stop (FINAL; STATUS.md "Coordinator ask 1 / 2", runs `*_p3`, job 55033059; VERDICT.md retires the 2 % cut): pursue
    while the best seed beats 2 SE and all 12 solvent controls, keep blocks that pass the same test left one out, and join blocks
    into one site when their supports overlap beyond hypergeometric chance + 2 SE. It stops by itself everywhere (blocks ->
    kept -> sites): 10 % fix 8 -> 8 -> 5 arms; 10 % null twin 14 -> 13 -> 5 moving arms (no discreteness: the leading blocks'
    2- vs 1-Gaussian gain -.0027 to -.0030 vs fix +.0016 to +.015); 5 % fix 7 -> 5 (4 arms + 1 weak block); qtoyflex 3 -> 1;
    5 % null 8 -> 0 and qtoynull 3 -> 0 (their seeds beat the controls, 5-37 nats vs <= 0.9, but the support-dof shell prior
    drains every grown block to zero and the leave-one-out check drops them). The 2 % cut's sixth 10 % candidate was a second
    seed on arm a5; without the cut two more blocks join a2, and the null's 8 extra blocks read no angle, pose, nuisance or
    quality. Blocks, not sites, carry states: a GMM on a site's union fails (full: fix a2 .003; tied: the null's 4-block a1
    site .011), and qtoyflex's one site has holo on its leading block (.928) and closed on its third (.900). Sites are masks;
    states are read per block. Open: which blocks of a site carry discrete states (the 2 % cut hid this on IgM; GAPw per p3
    block is igm-refs Part 4b, oracle session). Untested:
    the pursuit takes the best of ~60 seeds, so a best-of-N (Sidak) threshold of ~3.1 SE would trim the weak tail.
  - Final code (`x10fix_p2`) replicates the first run (purity .96-.97, supports 1,250-1,470 voxels, readouts within .01); 10 % null:
    8 blocks covering all five moving arms (purity .63-.92), continuous-angle R^2 .28-.31 as T1, BC .33.
  - Cost: 35-50 min per 100k particles on 16 cores, 31 GB; gate roots 5-7 min. The learned supports (`block%02d_support.mrc`) can
    serve as site masks now. Native port: per-block PCG supports + noise-calibrated group selection + the pursuit loop.
  - Lever 3 on these supports (igm-refs Part 3 above): the best chain measured, fully data-only (10 % product .338, 5 % .104).
  - Port (user: "Port sparse PPCA first", given to the sparse-ppca2 session): a native opt-in flex_pca path with a knob-free stop,
    reproducing the Python verdict, then (user: "After the port") the global-coordinate test at the r131 roots and the 692k root
    with |CTF|: does a site block isolate the closed arm better than the seam-free coordinate (.69-.72)?
    Design (approved 17:55; worktree `/Users/atifao/SIMPLE_sparse_ppca` at `e7f7efb80`; notes to `$H/handoff/sparse-ppca/PORT.md`):
    `flex_pca sparse_sites=yes` (default no; shared memory only). flex's prepare phase runs unchanged (sigma2, ml_reg, plane
    store, CTF or |CTF| per the stk flag). Three new modules in `src/main/flex/fit`:
    - `simple_flex_sites_ops`: the coarse-lattice operator, per-shell scale, noise and amplitude, the 6 rigid tangents as
      nuisance;
    - `simple_flex_sites_fit`: E-step, per-cell normal equations + ADMM group lasso, support-dof shell prior;
    - `simple_flex_sites_pursuit`: seeds and solvent controls, the knob-free stop, writers.
    PCG stays outside the pursuit (it would need one support per block); the site masks come out pcg_mskfile-ready. Outputs in
    `<run>/sites/`: site and block masks, per-block coordinates with the rigid tangents profiled out, nuisance coefficients; no
    state calls, no project writes. Working voxel ~lp/2.1 (the Python rule, not flex's lp/2.5). Validation at the Python lp per
    fixture (20 A IgM, 15 A qtoy), plus one IgM run at the tree's 15 A.
    Build (18:30): patch `$H/handoff/sparse-ppca/sparse_sites.patch` against `e7f7efb80` (8 files, two new modules
    `simple_flex_sites_fit` (operator + model) and `simple_flex_sites_pursuit`; the only added comments are the two !@descr
    headers); build `$M/sparse_port/b2` (job 55037253), fast gate 15/15, check_flex_dag and check_cli_chain pass. Native = Python on
    the gate roots: qtoyflex 2 blocks in 1 site, closed / holo AUROC .901 / .929 (Python .901 / .931), 118 s and 1.8 GB (Python
    ~5 min); qtoynull 0 sites. One seed constant (5 lattice voxels, ~2.3 lp) for both molecules. The per-round grown-block check
    stops qtoynull at round 1 and the 5 % null at round 2 (was 9 rounds), leaves qtoyflex's list unchanged.
    Native validation (5 of 6 fixtures, rule as written):
    - 10 % fix: 5/5 arms, 7 blocks in 5 sites, purity .95-.97, coverage .99-1.00, readouts within .006 of Python. 22.7 min and
      13.6 GB (Python 35-57 min, 31 GB).
    - 5 % fix: 5/5 arms (Python 4/5, a5 missed); R^2 .185-.304, LDA .102-.181. 26.6 min.
    - 5 % null: 5 sites, the five continuously moving arms (coverage .99-1.00, angle R^2 .135-.151, above T1's true-mask .12),
      with no discreteness (2- vs 1-Gaussian gain -.0028 to -.0033). The null twins' arms really swing, so Python's 0 sites
      there was a sensitivity miss: the criterion is "no discreteness", not "0 sites".
    - qtoyflex: closed / holo .901 / .929. qtoynull: 0 sites. 10 % null: running.
    Departure: under flex's noise model the solvent controls' held-out gains are negative (-1.6 to -2.9; Python +0.5 to +4.6),
    so near-dead blocks (0.3-2.2 nats, SE <= 0.4) pass "2 SE and every control" and widen arm masks.
    Fix, as delivered: a grown block must keep at least its seed's held-out gain (strict; the 2-SE tolerance was asked for but
    not built). Real blocks end at 1.04-16x their seed and weak ones at 0.07-0.67x. The floor removes the near-dead blocks
    (four on the 5 % null, one on the 5 % fix) and no real block on any fixture. A block that is strong when added and redundant
    later stays inside its arm's site, and GAPw rejects such blocks (Part 4b).
    COMMITTED on `hetero-tree-wip` as `b3ced0021` (build b3, fast gate 15/15; `$H/handoff/sparse-ppca/PORT.md`):
    - every readout within .01 of the prototype;
    - 10 % null: 5 arm sites, no discreteness;
    - 18-40 min and 13.6 GB per 100k IgM particles (prototype 46-125 min, 30.7 GB); request 24 GB per 100k at box 70.
    Real-data test running on b3: r131 roots on both seeds (55039094_1-2) and the 692k root (55039095_3, 64 GB). They report
    closed AUROC / d' and the closed share of the top-N tail per block, and list every block the strict floor drops.
    - r131 RESULT (b3, lp 15, box 36 at 6.5 A, |CTF| from flex's planes; both seeds share the root; `native/rr131s2`): no site
      block isolates the closed arm better than the seam-free coordinate.
      - The pursuit accepts 18 blocks and stops at round 19; leave-one-out keeps all 18. The molecule is compact, so all 18
        overlap and form one site.
      - Cross-half closed AUROC (d') per block: best .659 (.58), then .597 and .574; the other 15 .50-.55.
      - Closed share of each block's top-3,075 tail: .072-.111 (base .069).
      - The whole site's span (36 coordinates): AUROC .711 (d' .78), tail share .098. That reproduces seam-free's .69-.72 and
        does not beat it.
      - So RNase P's closed signal is spread over the compact molecule, not held by a localised site, unlike the IgM arms.
    - TRAP (confirmed, PORT.md; diagnostic run `native/rr131diag`): on r131 any local mode captures ~100 held-out nats of
      reproducible non-site variance (noise-model or per-particle misfit), wherever it sits.
      - mskdiam 225 A in a 234 A box leaves room for ONE solvent control (at the rim, 112 A from the centre), not 12. It gains
        109.5 nats (SE 9.9).
      - The 37 on-molecule seeds gain 72-275 (median 140) with no radial trend (r -.14).
      - By round 19 the unused seeds' median (69) falls below the control (90), and the pursuit stops there.
      - The 692k root shows the same floor per particle.
      So on real data the stop is set by that misfit floor, and blocks near it may be misfit, which fits no block isolating the
      closed arm.
      - The strict seed-gain floor dropped nothing on r131: all 18 blocks kept 1.5-4x their seed's gain.
      - Fixes to weigh before using sparse sites on real data: more controls than a tight box allows, or modelling the
        misfit (e.g. per-particle band-wise scale).
    - 692k RESULT (2 h 10 min, 15.6 GB, one node; 8,720 closed reference rows):
      - Pursuit: 16 blocks, one site. It stops at round 17 when the best seed (1,830 nats) falls just under the single solvent
        control (1,834). The misfit floor per particle is about r131's.
      - Closed vs rest (AUROC / d'): block 3 alone .802 / 1.17, at 67 A from the centre; the span of all 32 coordinates .849 / 1.43;
        seam-free z_MAP .798 / 1.13; z_ML .667 / .45. So one 2-coordinate block matches the seam-free coordinate, and the span
        beats it modestly.
      - Purity stays low. The top of block 3's leading axis is 5.3 % closed (4.2x the base rate), while the tree's best leaf
        `r1111221_s2` is 9.3 % closed.
      - The closed leaf pair (`pair_rank.py`, the seam-free test's own pair scorer on its 112 pairs): the span separates it more
        (d' .48 vs .25) but separates every pair more (median .61 vs .23), so it ranks 75 / 112 (seam-free 46). No single block
        singles it out (best: block 14, 18th).
      - The strict seed-gain floor dropped no block at either root (692k: every block kept 1.87-4.19x its seed's gain).
      - Reading: on RNase P the sparse sites give a slightly better closed-vs-rest coordinate at 692k and none at r131. They do
        not single out the tree's closed leaves, and their tails are less pure than the tree's best closed leaf.
  - Side finding: a converged zero-penalty whole-molecule fit reads a1 / a3 (span LDA NMI .18 / .19) where flex's whole-molecule
    basis reads <= .02; flex's basis stays near its data-free start (2 + 1 iterations).
- `spectrum-stop`: a stopping rule from each node's flex spectrum (directions above a noise edge set by the node's N, after the
  spiked-covariance phase transition of Baik, Ben Arous & Peche 2005; Katsevich, Katsevich & Singer 2015 in cryo-EM), and a
  global bound of r + 1 identifiable discrete states for r significant directions. VERDICT (`$H/handoff/spectrum-stop/VERDICT.md`,
  433 nodes of 44 trees): NO. Nulls carry above-edge nuisance directions (qtoynull 1.75 vs edge 0.77: pose and contrast; the IgM
  fixture and null twin have identical spectra led by shift-error tangents, arm spikes <= 0.08), and the RNase P closed spike is
  far below any edge (0.014 at the r131 root, 0.01-0.12 down the 692k closed lineage), so the rule would fail the hard gate; a
  null-calibrated edge stops the nulls but loses qtoy2's close pairs, the qtoyflex root and r131. At 692k the stops track N, not
  exhausted variance. The global r + 1 is only a lower bound and fails where states exist at r = 0 (r131; IgM, 32 combinations).
  The empirical rare-state limit (n/sqrt(N) >= 16-17) IS the BBP law, with p_eff 1,000-2,600: for the paper, not a rule.
  Site-focused latents raise the arm spikes 15-140x, but the null's continuous arms spike too.
- `flex-precision` (verify-flex-pca-s-embedd session; runs `$M/flex_precision/`; VERDICT and both patch variants to
  `$H/handoff/flex-precision/`): CONFIRMED in source and numerically. The embedding stores the MAP sampling precision
  Q = A G^+ A (rebuilt from the captured per-particle Gram blocks: equal to Q to 4e-15, 99 % off A; prior-whitened eigenvalues
  never below 4 on 60,000 particles), while `calibrate_noise_scale` / `deconvolve_latent` read it as A and `mest_charts.py` as G.
  43 % of the per-particle directions are weak (data / prior < 1). On its test node flex calibrates the noise scale to 0.460
  (fixed reading 0.507, truth 0.634), so the deconvolution under-scales the real conformer separation by 26-42 % on the strong
  axes and turns noise into apparent signal on the weak ones (2-9x). Consumers: state seeding, site_rank starts, the charts code.
  The patch touches the import line and the three `map_sampling_precision` calls in `simple_flex_pca_embed.f90`; it applies
  stacked on the pose-profile patch (`flex_precision_posterior_on_pose_profile.patch`), not in the other order with plain
  `git apply`. VERDICT: storing A instead of Q fixes the strong axes (z1-z3, 84 % of the truth variance: deconvolved variance
  26-42 % low -> 10-20 % high) and the per-particle posterior means (error .385 -> .299 of the truth variance, -22 %), but not the
  weak-axis inflation, which comes from the half-split noise scale (.507 vs truth .634; likely error both halves share).
  Gated (VERDICT.md "Gate"; build `$M/flex_precision/build_fp` = `e7f7efb80` + the stacked patch; arms `fprec` vs pp_on / `fpbase`):
  fast gate 15/15; nulls one leaf; the r131 corner kept on both seeds (1.01-1.02 / .57; pp_on .52-.54 against the .6 limit); qtoy2
  NMI .760 / .799 vs .811 / .813, with closed merged into combo 22 on seed 1 (the root split into 4 groups instead of 5; pp_on kept
  closed on both). At the root's first flex call the fix moves labels against truth from NMI .456 to .512 and the class separation
  the coordinates keep from .48 to .70; the noise-scale gap is not pose error (1.25 with the profile on or off). A consumer the
  first audit missed: the root's first flex call's states build fflex's focus mask, so the fix also changes fflex's basis
  (eigenvalues up to 12 %). Third qtoy2 seed pair: the fix merged closed again (NMI .753 vs .810), so 2 of 3 seeds vs 0 of 3
  without it. HOLD, not committed.
- `hier-ppca`: the tree as a hierarchical PPCA (after Bishop & Tipping 1998): every node keeps its PPCA, children are carved in
  the parent's latent and refined, a view-matched map test selects, a continuum stays one leaf with its own coordinate.
  VERDICT (`$H/handoff/hier-ppca/STATUS.md`, scripts beside it, runs `$L/hppca/<key>`):
  - **A cross-half map test of hard-assigned groups is valid only on data no assignment step used.** In the band the assignment
    scored, the 5 groups of a one-conformer null (qtoynull) are reproducibly distinct at z 53-56, with cross-half (swap) scoring
    and 300 view x defocus cells, and with no view, defocus, SNR, B or pose-error difference in truth: each half's groups were
    selected by references carrying the other half's noise. Band split (assign below lp L, test only above it): 5 -> 1 at L 10
    (test 9.4-5.6 A, max z 1.9) and at L 15 (14.6-5.6 A, 6 shuffles, max z 1.97).
  - The production split test does not meet this (checked in the source): `state_split` runs its cross-half EM below EM_LIM 20 A
    and tests from 4 shells above it (the tree passes a focus mask) to lp 10 A, starting from stage-B labels assigned at LPB 10 A, so its test band sits inside
    the band that chose the labels. The gate nulls pass through `locality_guard` + `effect_floor`; the IgM nulls do not
    (10 % 50/50, 5 % 64/36).
  - With the band split the select step keeps real states apart, but the model does not beat the tree: qtoy2 root 8 leaves,
    NMI .697 (stage B .665, lever 4 .705, tree mm2 .780), closed F1 .916, holo .818, the pairs 1/2 and 43/49 merged; qtoyflex holo
    isolated (F1 .689, tree .62) but closed lost (.22, tree .39-.47).
  - Continuum vs discrete is not decidable from first-moment maps at these SNRs: in the assignment band the sorting ramp
    dominates (the IgM fix and null arms read the same), in the clean band 8 bins have no power and 4 bins break at the group
    boundary for a continuum too; single linkage keeps the qtoyflex arm together but bridges qtoy2 classes (NMI .54).
  - IgM 10 %: no call; the arm information lies below 15 A, so a band split leaves the test nothing.
  - Nested coordinates (qtoyflex, s1 R^2): root latent .187, level-1 leaf .289, leaf + local .314, + root latent .339.
  - Traps: per-shell amplitude matching is degenerate when one map has no band signal; groups under ~2k particles have no
    clean-band signal (flag half-map correlation < .05 as untestable); 2 shuffles per pair give 1 df (use >= 6).
  - Next for continuum vs discrete: gold-standard classification (independent half refinements, half-specific flex fits, groups
    matched across halves), where the test can use the full band.
- `fission` (the f2r session; `$H/handoff/fission/STATUS.md`, runs `$L/fission/`): data fission (Leiner et al. 2023; Neufeld et
  al. 2024) as the alternative that keeps the low band: X1 = X + tau Z labels, X2 = X - Z/tau tests, Z white in real space at
  each image's noise sd (stacks `$L/fission/stacks/<toy>_tau<t>/X{1,2}/all.mrcs`). Validated: no cross-copy noise correlation
  (-.0004 / +.0002); at tau 1 each copy has twice the noise variance, at tau 0.6 X1 keeps more (corr with X .86) and X2 pays (.51).
  qtoynull tau 1: the over-split K 4 on X1 merges 4 -> 1 on X2 in lever 4's full band (max pair z 0.46, z_crit 3.38). qtoy2 tau 1: the X2 test at
  39-10 A keeps every real state (pairs of different composition z >= 6.9; a clean 28|28 duplicate merges at z 1.5); labelling on
  X1 costs .07-.09 class accuracy (.703 vs lever 4's .785 on X; closed F1 .83 vs .95). A holo view split (a tilt 105-150 cone)
  survives at z 8.9: the 120-cell shuffle cannot mimic a split finer than its cells (480 / 1920-cell rerun: job 55031241). The
  9.4-5.6 A band has no power on small or low-SNR states (it folds closed into a junk lump).
  Power (hier-ppca, qtoyflex truth groups, z): closed vs arm 17-22 unprotected, band split 0.5-3.4 (lost: the difference lives
  below 15 A), fission tau 1 9.3-14.8, tau 0.6 4.6-8.9; null .98 / .04 / .75 / .30. X2's agreement follows S/(S + cN) with
  c = 1 + 1/tau^2. Recommendation: fission at tau 0.6 for the select step (X1 keeps more SNR: qtoy2 NMI .633 after round 1 vs .570
  after 6 rounds at tau 1), tested at 39-10.2 A on X2, groups under ~2k untestable; the band split only where the noise model
  cannot be trusted. Real data need an exact per-particle noise model (coloured, per micrograph; r131's registered stacks are
  non-stationary: extremes-repro).
  tau 0.6 (55030590-1): the null still passes (X2 4 -> 1, z <= 2.7; X1 4 -> 4 at z 5-11), and labelling on X1 recovers most of the
  cost: qtoy2 class accuracy .751, closed F1 .904 after 6 rounds (tau 1 .695 / .803; on X .772 / .947); X2 at 39-10 A keeps all
  10 states (the holo pair at z 4.9). Finer nulls (tau 1, 480 and 1,920 view x defocus cells) keep the holo pair (z 8.4 / 8.8), so
  it is not only a coarse-strata artefact; a control job (55032828: true-holo restriction, view-cone split of true holo, non-holo
  remainders) separates the pair's admixture from its views.
  VERDICT (`$H/handoff/fission/VERDICT.md`): fission removes the self-selection. Labels from X1, maps from X2: the qtoynull root
  merges back to one state in both bands, with either seed, at both tau; qtoy2 keeps every state of different composition
  (z >= 6.9 at tau 1, >= 4.1 at tau 0.6) and merges a clean duplicate. Labelling cost vs the same passes on X: .08 class accuracy
  at tau 1, .02 at tau 0.6 (.751 vs .772, closed F1 .90 vs .95). The one kept same-class pair is a purity split, not a view split:
  its pure-holo cores are equivalent (z 2.1-2.5) while its non-holo remainders (18 % vs 10 %) differ (z 6.0), which no map-level
  test can merge. On real data Z must be shaped by the per-group noise spectrum (the committed sigma2 groups or flex's noise model).
- `fission-split` (fission-split2 session; `$H/handoff/fission-split/STATUS.md`; worktree `/Users/atifao/SIMPLE_fission_split` at
  `cc6d9cbb6`, prototype env vars only):
  - Setup: the whole tree runs on an X1 project, so poses, labels, the site layer and state_split's EM never see X2. Only the
    split decision reads X2 (`SIMPLE_SPLIT_FISSION=<X2 project>`), and `SIMPLE_SPLIT_GUARDS=no|locality` switches the guards.
    tau 0.6. The copies are validated by the fake-map null (qtoynull .010, olenull .011, micnull .015, r131 .025).
  - (a) Today (mm2, pre-guard log): all three gate nulls pass the split test (qtoynull z 4.68, olenull 9.20, micnull 75.1),
    and locality_guard merges each.
  - (b) Guards off, fission test, same labels read on X1 vs X2:
    - qtoynull MERGED (X2 z -0.19 vs X1 4.10): self-selection.
    - micnull MERGED (X2 z 0.41 vs X1 61.0): self-selection with a weak pose tilt (eta^2 .023), not the pose-scatter split
      expected.
    - olenull KEPT (X2 z 7.78, 2 leaves; X1 15.2): a real contrast / brightness difference. The cross-half EM's small group is
      brighter in truth (snrq 1.19 vs 1.09); B, defocus and pose are equal. The locality guard merges it (L/L_ref .228).
    - The effect floor is redundant under fission on all three.
  - So the guards were doing two jobs. Fission replaces the self-selection one; the locality guard keeps the "states differ
    locally" one.
  - (c) Fission + the locality guard only, build fsp_b2 (patch v3, fast gate 15/15):
    - qtoy2 (55038253): the root keeps 5 groups (weakest pair z 22.5). At depth 1, r2 (z 15.9) and r4 (11.5) are kept and r3
      (1.2) is one leaf.
    - r131 FAILS THE HARD GATE at tau 0.6 on both seeds (55038254-5): the root merges to one leaf of 44,334, so the corner is
      never reached.
    - Diagnosis at the root's own labels (z_crit 3.71):

      | seed | X2 test (the tree) | X1-map replay | pp_on on X |
      | --- | --- | --- | --- |
      | 1 | z 1.97 | z 10.6 | z 15.8 |
      | 2 | z 2.81 | z 13.8 | z 18.5 |

      About 50 % (seed 1) and 70 % (seed 2) of the X1 between-group signal is real; the rest is imprint. The real part is too
      weak for X2's ~26 % of the information at tau 0.6, so the merge is lost power.
    - TRAP: re-reading groups selected on X on X2 is not honest, because X2 = S + N - Z/tau shares N with X. Only selection on
      X1 makes an X2 read honest.
    - tau 1 (55038926-7; the extremes session's tau 1 copies, which pass the fake-map null; projected z ~4.4 / 6.0): r131 STILL
      ends as one leaf of 44,334 on both seeds, so the hard gate fails at tau 1 too. qtoynull and micnull at tau 1 stay one leaf;
      olenull at tau 1 is running.
    - Reading: an honest split test cannot confirm r131's root split at 44k particles at either tau. That agrees with
      extremes-repro: the corner passes partly on imprint. By the coordinator's rule the session stops here (no other tau or
      knobs) and writes the verdict.
  - VERDICT (`$H/handoff/fission-split/VERDICT.md`; patch v3 against `cc6d9cbb6`, not committed):
    - Fission + the locality guard passes every gate null at tau 0.6 and tau 1. qtoynull z -0.19 / 0.42 and micnull 0.41 / 0.80;
      olenull keeps a real contrast split at tau 0.6 that the locality guard merges, and is one leaf by fission alone at tau 1.
      The effect floor never acts.
    - qtoy2 at tau 0.6: 8 leaves, NMI .717 (pp_on .811); the 1/2 and 22 pairs merge.
    - r131 FAILS the hard gate at both tau: X2 z 2.81 / 1.97 at tau 0.6, 0.88 / 1.67 at tau 1 (z_crit 3.71).
    - The r131 root is not the closed split in any tree (closed enrichment of root states x0.87-1.13). It is a near-tie: root NMI
      between seeds .43 on X, .04 under fission.
    - The noisier tau 1 labels change the partition itself (the real share of X1's signal falls from 50-70 % to 18-22 %).
    - Honest X2 rebuilds of the root maps read .83-.87 closed: no corner. The corner sits under a weak root near-tie that no
      honest test confirms at 44k particles.
    - Untested routes that keep full-data labels (band split, gold-standard halves) are listed in the VERDICT's last section.
- `extremes-repro` (runs `$L/extremes/`; notes to `$H/handoff/extremes/`): whether the extreme maps along a latent coordinate
  reproduce across fission copies. qtoyflex (first numbers): self-selected tail maps "reproduce" across halves (.13-.51 at
  39-10 A), but so do tails picked along a fake map no particle carries (.23-.45); honest tails (select on X1, rebuild from X2)
  reproduce at .00-.09, and fission removes 98-99.7 % of the self-selection imprint. Self-selected extreme maps read 2-4x further
  out along the true continuum than their particles are; honest ones track the true positions, at a cost (X1's half SNR: 5 % tails
  .08-.11 of the arm range vs .23-.28 when X selects). Cross-half reproducibility of self-selected extremes is no evidence.
  Traps: r131's registered stacks have non-stationary noise (corners .36-.66 of the inside low/mid-frequency noise; a band-wise
  radial model fixes it); the map-level even/odd cross-copy test cannot validate fission on real data (only the fake-map null
  can); the CTF-flip trap (section 9) was found here.
  VERDICT (`$H/handoff/extremes/VERDICT.md`, METHOD.md, RUNS.md, tables/<node>/): self-selection is large for tails and leaves
  alike: at the same particles, own-data maps place extremes 1.3-3.1x further along the true motion (IgM arms 1.74-2.5x even for
  90-98 %-pure tails; leaves 1.75-1.8x), and fake maps "reproduce" at .16-.45 own-data, 0 honestly. Fission removes it (fake-map
  leak <= 1.8 % qtoyflex, <= 3.5 % r131 with the band-wise noise model); honest maps read the truth on one scale (IgM a1 beta /
  shift 345-360 at every tail size vs 613-793 own-data). r131: the corner leaf r2_s2 reads 1.03 / .54 own-data, .93 / .59 honestly;
  no leaf passes the hard gate honestly, and only small latent tails reach its threshold honestly (seam-free pc2 bottom 5 %
  .966 / .574, site-2 focus top 5 % .953 / .579, ~2.2k particles). qtoyflex: the tree beats latent tails at matched N (closed
  .197 vs <= .13). IgM 10 %: focused coordinates on a data-only site mask find what the tree misses (a1 tails 90-98 % pure vs 62 %
  for the leaves; an a2 tail 87 % in its -30 deg state with no a2 leaf; null extremes +-12-15 deg vs +-1.7). So neither the global
  coordinate nor the node latent identifies extremes better than the leaves; focused site coordinates do for small parts, and
  every own-data readout, the r131 gate included, needs a fission re-read. Caveat: the honest leaves are a re-assignment among
  the tree's leaf references on X1 (agreement with the delivered labels .50 r131, .61 qtoyflex, .62 IgM fix), not a tree run on X1.

### 7.9 A cohesive locally linear / mixture model: ideas from structure from motion and related fields (2026-11-07)

User: "I would really like a cohesive locally linear or mixture ppca ... look into more literature on structure from motion or
fields with similar problems". Literature scan by a subagent (42 searches; citations checked where marked in its report):
- Non-rigid structure from motion:
  - Torresani, Hertzmann & Bregler 2008 (EM-PPCA at known motion = our model).
  - Procrustean normal distribution (Lee et al. 2013; mixtures: Cho, Lee & Oh 2016): deformation in a rigid-aligned frame, as
    our pose profile does.
  - Russell, Fayad & Agapito 2011: overlapping local models, fitted as a multi-label energy whose overlaps must agree.
  - Yan & Pollefeys 2008: linked parts' motion subspaces intersect, so kinematic chains come from subspace ranks.
  - Unions of subspaces (Zhu et al. 2014; Kumar et al. 2017-2019; Agudo 2022, discrete classes plus deformation).
- Coordinating local linear models (Roweis, Saul & Hinton 2001; Brand 2002; Teh & Roweis 2002; LTSA 2004): seam-free by joint
  alignment. At d' 1-2, though, neighbouring discrete states share responsibility and get glued into one continuum.
- Mixtures from projected data:
  - Common-loading factor-analyser mixtures (Baek, McLachlan & Flack 2010) and extreme deconvolution (XD; Bovy, Hogg & Roweis
    2011).
  - Compressive classification (Reboredo et al. 2013; Yu, Sapiro & Mallat 2012): a minority is separable only along P_i Delta
    outside the projected continuum covariance. A 5 % state adds ~.05 |Delta|^2 of variance, so variance ranking (BBP) can
    bury a state that likelihood still separates.
- Single-cell trajectory inference (PAGA, Slingshot, DDRTree, Mapper, cluster trees): connectivity is decided on per-sample
  neighbourhoods, which call discrete states connected at d' 1-2.
- Model order reduction: quadratic manifolds (Barnett & Farhat 2022), x ~ x0 + Vq + W(q kron q). They stay additive in density
  (unlike the falsified deformation field) and curve a low-dimensional coordinate.
- Cryo-EM:
  - RECOVAR (Gilles & Singer 2025): fixed-pose PPCA + deconvolved density + paths; lambda by L-curve; <= 4-D.
  - CryoBench 2024: RECOVAR is the best fixed-pose method on IgG-1D, Spike-MD and Ribosembly.
  - Hydra (K neural fields x d latent); 3DVA; Andén & Singer 2018.

Against our record:
- The field's cohesive answer is a shared basis with a mixture in latent space (common-loading mixtures / XD). flex_pca already
  does this for the whole molecule (deconvolved latent -> XD labels), and the flex-precision bug sat in its input.
- Our verdicts say the whole-molecule version cannot see the states:
  - the chart XD gap is blind below d' ~2.5;
  - latent density basins have no power;
  - the global coordinate isolates no minority.
- What breaks the SNR barrier here is locality: sparse PPCA's BBP-local supports.

Proposal: one cohesive generative model per node, y_i = C_i P_i (mu + sum_s B_s z_is + T n_i) + noise. The site bases B_s are
sparse PPCA's learned supports, the rigid tangents T are profiled, and one joint latent density over the site coordinates is fitted
with the posterior precision A. The tree stays as the recursion that refits poses per subset. Pieces, cheapest first:
1. Soft XD over a site's joint block coordinates as the seed. It can supply the split direction the saturation test needs for a
   minority next to a continuum (qtoyflex's closed state sits on its site's third block; round 1 missed it). Test: qtoyflex closed
   recall and purity, then Part 5's call; IgM 10 % fix / null as the regression.
2. Cross-site coupling (Yan & Pollefeys): the rank of the deconvolved cross-site latent covariance separates factorial
   (independent) sites from coupled ones, and a coupled direction is where a global coordinate comes from. Test: IgM (independent
   by construction) and a coupled synthetic; r131 after the port.
3. A quadratic manifold per site for curved swings: one q per site, with W (q kron q) as extra basis volumes in flex's normal
   equations and q on a grid. Test: qtoyflex arm, R^2 against true s1, closed off the curve.
4. Per-site decisions by the fission saturation test (igm-refs Part 5), not by latent density gaps.
5. Minority directions ranked by separability, not variance (compressive classification): the RNase P closed spike sits far below
   every BBP edge.

User, same evening: "a cohesive model over the entire tree". So the site basis is shared by the WHOLE tree, not refitted per node,
and the mixture over it is tree-structured. Tests 1-3, re-aimed and started as sessions:
- `tree-latent`: the tree built in latent space over one shared site basis, splits accepted by the saturation test; scored against
  truth, the hetero tree and lever 3. Notes to `$H/handoff/tree-latent/`.
- `site-coupling`: factorial (a product of per-site trees) vs one joint tree, from the cross-site coupling. Notes to
  `$H/handoff/site-coupling/`.
- `quad-manifold`: one curved coordinate per site. Notes to `$H/handoff/quad-manifold/`.

tree-latent interim (19:20). The tree grows on fission copy y1, since membership is a selection too, and splits are accepted by
the saturation test.
- qtoyflex: 2 leaves. The root split pulls closed and holo away from the arm continuum (DISCRETE, z 5.2). Closed's best node
  is r1, a depth-1 leaf shared with holo: purity .21 / .24 and recall .34 / .54 (y1 / y).
- qtoyflex, closed vs holo below r1: the split was proposed correctly but gets no call at 5.8k particles.
- qtoyflex scores: kind NMI .325 on y vs the gate trees' .255 / .265; the arm stays one leaf. So the latent XD split supplies
  the minority direction that round-1 labels lacked.
- IgM 5 % fix: one leaf (joint tree and forest alike). On y1 alone the deconvolution cannot resolve d' ~1 arms: the same
  fission power cost as Part 6.
- Running: IgM 10 % fix / null and qtoy2, the forests, and a v3 with proposals and labels on y and membership on y1.

quad-manifold VERDICT, half yes (`$H/handoff/quad-manifold/VERDICT.md`, tables.md; runs `$M/quad_manifold/runs`):
- A curved swing gets one good coordinate. IgM 10 % null a1: the 1-D curve with a learned density ties the 2-mode block
  (+12 +- 62 held-out nats). It beats a straight line with the same density freedom by +3,461 +- 111, reads the angle at R^2
  .342 (block .287, line .304), and stays monotone where the line folds.
- States do not sit off the curve. With a learned density the likelihood bends the curve to take the strongest state as a mode.
  qtoyflex holo becomes the curve's far end and closed merges with the arm's end (off-curve AUROC holo .40, closed .53; s1 R^2
  .295 on the curve vs .347 line, .339 block).
- Curvature does not change the discreteness call. The saturation test on q-quantiles gives fix DISCRETE (.70 / .68 curve,
  .67 / .65 line), null not discrete (2.95 / 2.59, 2.35 / 1.55) and qtoyflex not discrete.
- The curve-minus-line held-out gain is a curvature readout: +3,461 on the swing, -167 on the 2-state fix. Three non-collinear
  states would also read "curved" (untested).
- With per-particle b_i, G_i, the 2-atom vs segment readout separates both IgM twins (fix +561, null -3,385), unlike Part 4.
  Plausible, but no homoscedastic control was run.
- Port traps: every site model needs a free reference offset c0 with its own prior (the qtoy gate roots' mean misfits a site by
  4.1-4.8k nats), and a 1-D curve fitted to 2-D residue folds and shows false modes.
- 2-D (final; VERDICT.md item 7, STATUS.md addendum): 2-D q follows the 1-D pattern, so curvature pays only on the curved
  swing.
  - IgM null a1: +3,341 +- 96 over the block, of which curvature (vs a 2-D linear model with the same learned density) is
    +2,846 +- 92.
  - IgM fix a1: +2,059 +- 77, all of it density (curvature -9 +- 55).
  - qtoyflex: ~+4,000, about half each (log plateau only; its 2-D runs died of lattice growth).
  - The second dimension adds on every structured site, so none of them is 1-D.

site-coupling VERDICT (`$H/handoff/site-coupling/VERDICT.md`): factorial across sites, joint within a site.
- IgM product-form subset: per-site K comes out 2, 2, 2, 3, 3 (the truth's 72 cells). The five-arm product beats a free
  72-cell table by 45 held-out nats (z -5.3) and a joint mixture by 2,057.
- Full sets: the product beats the joint mixture by 2,030 nats (fix twin) and 1,651 (null twin).
- Planted couplings are found: a1-a3 at z 15.6, coupled null angles at z 9.6. No false calls on 65 independent pairs.
- qtoyflex's one site is coupled internally.
- A covariance screen alone false-calls 6 / 90 null pairs, so it can shortlist but not decide.
- TRAP: integrating nuisances under their fitted prior invents couplings; profile them.
- Model shape: one tree per site (a particle's leaf = the product of its site leaves), a joint tree over each site's blocks; two
  sites merge only when their free table beats the product by 2 SE held out. This is also why the joint hetero tree needs an
  equivalence merge-up: it duplicates a site's states under other branches.
- Sent to tree-latent (the port note builds on the forest) and site-chain (a factorial output-mode proposal after its tests).
- Final (all jobs done):
  - The null twin's 10 pairs all favour the product (z -1.1 to -7.5); 0 false calls on 70 independent pairs.
  - CAVEAT: refits of the site mixtures from another start keep the decision but move the gain by up to 5 paired SE, so fit
    site mixtures with restarts before the merge test.
  - The covariance screen catches 2 of 4 same-arm IgM block pairs.

tree-latent, qtoy2 v2: the root split is DISCRETE (z 28.8, deconvolved d' 9.9), but the numpy code took 51 min for qtoy2's
14-coordinate site; a compiled version is rerunning.
tree-latent, IgM 10 % fix v2 (joint tree on y1, frozen poses):
- 6 leaves; five-arm product .223 on y (.135 on y1), against the v2 forest's .208, lever 3's .338 and the hetero tree's .042.
- Per arm: a1 .234, a3 .241, a4 .023, a5 .039.
- So at half the information it beats the tree by far but stays below the masked-round chain. v3 (proposals on y) is pending.

tree-latent VERDICT (`$H/handoff/tree-latent/VERDICT.md`). The candidate cohesive model is a tree grown in latent space over
ONE tree-wide sparse site basis: deconvolved node Gaussians (a Bishop-Tipping soft hierarchy with XD, posterior precision A),
growth by the saturation test, frozen poses. In v3, proposals and test labels are on y; memberships and child weights are on y1.
- qtoyflex: 3 leaves (arm kept whole, holo, closed). Closed is its own depth-2 leaf (F1 .558, purity .63, recall .50); holo
  .713; kind NMI .399. Gate trees: closed .39 / .47, holo .62, kind .255 / .265.
- qtoy2: 8-9 leaves, NMI .752-.764, closed .84-.90, holo .95-.96. It loses the conformer 1 | 2 pair. pp_on .811.
- IgM 10 % fix: five-arm product, joint .260 / forest .258 (tree .042, lever 3 .338). The 2-state arms reach lever-3 quality.
- IgM 5 % fix: 1 leaf (root z_mid 2.2-2.8: fission's power cost at d' ~1).
- Nulls: 1 leaf everywhere, including the IgM 10 % twin that the hetero tree splits. No null accept in any variant.
- Forest = joint. Balanced K = 2 splits of a 3-state arm straddle and read continuous; lopsided ones pass.
- Next (approved): K = 3 proposals with a pairwise saturation test in forest form (for a2 / a5 and the qtoy2 1 | 2 pair).
- K = 3 (interim, STATUS addendum):
  - At IgM 10 % the root recovers the 3-state arms' proportions: a2 has both boundaries DISCRETE; a5's 20 % state splits off,
    while its two 40 % states read continuous.
  - qtoyflex keeps the arm whole and splits closed | holo. IgM 5 % stays one leaf.
  - Null margins get thin under K = 3 (the closest null side +.043 at 2 SE). A Sidak correction for the side tests per node
    restores +.083 and so far costs nothing on a2 / a5; it drops one spurious a1 cut.
- PORT question (user: "port the shared tree latent if it's useful yet"). Useful on the synthetic fixtures, but untested on RNase
  P, where the r131 corner is the hard gate an honest test already failed (fission-split). Before any Fortran, the prototype
  runs on the native sparse fits of the r131 root and the 692k root, at frozen poses, and reports:
  - the closed share per leaf, and closed's own node;
  - the hard-gate readout on honest y2 maps;
  - each accepted split's correlation with per-particle nuisances, since real-data misfit (PORT.md trap) survives fission.
  A port design note comes in parallel (inside hetero_tree with per-subset pose refinement; cost at 692k).

Prior art for one shared tree model (user: "our current approach lacks elegance ... surely other people have done this before";
searched 2026-11-07):
- Cryo-EM: no single tree-structured generative model with per-node pose refinement found. The cohesive models are flat:
  - cryoDRGN: one decoder and a continuous latent; hierarchical classes come post hoc;
  - Hydra (Levy et al., NeurIPS 2024): a flat mixture of K neural fields, with poses, conformations and classes jointly;
  - hyper-molecules (Lederman & Singer 2017; Lederman, Andén & Singer 2019): tree topologies are only mentioned;
  - RECOVAR / 3DVA: one linear model at fixed poses.
- Elsewhere:
  - Bishop & Tipping 1998: a hierarchical mixture of latent variable models, soft, top-down EM.
  - Tree-structured stick breaking (Adams, Ghahramani & Jordan 2010): data at internal nodes, node parameters diffuse from the
    parent.
  - Phylogenetic factor analysis (Tolkoff et al. 2017; Hassler et al. 2022): latent factors diffuse along a tree by Brownian
    motion; one Gaussian model, internal nodes inferred.
  - TreeVAE (Manduchi et al. 2023): grows a tree inside one generative model.
  - Partial common principal components (Flury 1987) and shared-subspace models (Franks & Hoff): shared axes plus group-specific
    ones.
  - Latent tree graphical models (Choi et al. 2011).
  - The hierarchical Dirichlet process (Teh et al. 2006: shared atoms, group weights; from memory, not checked).
- Our version: y_i = C_i P_i(theta_i) [mu + B (m_node + z_i)] + noise.
  - B: one tree-wide sparse site basis; it grows, as in partial CPC, when a subtree finds a site.
  - m_node: node means that diffuse down the tree; tree_wiener is the posterior of that prior.
  - z_i: within-node coordinates. Continua live here, and a particle may stay at an internal node, as in TSSB.
  - theta_i: poses refined per subtree, which is block-coordinate ascent of the same likelihood.
  - Growth: by the fission saturation test.
- Pieces of the current tree it absorbs:
  - per-node flex fits -> B and m;
  - the site layer -> B's sites;
  - the guards -> the saturation test;
  - tree_wiener -> the diffusion posterior;
  - the equivalence merge-up -> a factorial forest. A joint tree duplicates one site's states under other branches;
    site-coupling decides between the two.

Not to retry:
- map-level soft EM (noise-locks);
- post-hoc chart coordination (seams);
- latent density basins / kNN connectivity (no power at d' 1-2);
- a deformation field replacing the basis (falsified; composition);
- VB/ARD mixtures.

### 7.10 Overnight plan, 2026-11-07 (user: "please work autonomously on the stuff overnight")

REFINED PLAN (2026-11-08 ~00:30; each thread has "if -> then" rules the coordinator applies without asking):
1. RNase P noise model (rnasep-misfit; the critical path). It tests a per-particle band-wise noise variance (empirical Bayes),
   an off-mask estimate and per-micrograph sigma2 groups. Each is scored on the r131 solvent floor (109.5 nats now), the
   y1-y2 correlation by noise decile, closed AUROC kept, and IgM as the regression.
   - If one wins on r131: run 692k.
   - Then tree-latent reruns its decisive RNase P test, with the noise model in the whitening AND the fission draw, and
     sparse-ppca2 reruns r131 / 692k.
   - If the shared model then finds the closed lineage honestly, the forest port starts (`handoff/tree-latent/PORT_DESIGN.md`).
   - 01:10: the band-wise noise variance removes most of the floor. Runs `$M/rnasep_misfit/runs/floor_r131_*`; each number is
     a block's held-out gain on one of 24 solvent control balls.
     - The 12 balls outside the support sphere (box corners and edges): base 34-133 nats. nvarB (per-particle band-wise noise
       variance, empirical Bayes) gives -5 to +8, and nvarBC -7 to +2.
     - The 12 in-sphere solvent balls (the native run's rim control): base 175-307. nvarB gives 7-50, nvarBC -3 to 35 and
       insideB 4-44.
     - A scalar per-particle variance (nvar1) only reaches 129-258 / 22-98, and a per-particle band scale (bscale) changes
       nothing.
     - So the floor is the noise level, band by band. Up to 35-50 nats remain on the in-sphere solvent.
     - rnasep-misfit already runs the shared tree with nvarB, base vs nvarB, on r131 and a 692k 50k subsample (jobs 55052871,
       55054551). `mf_tree.py` is tl_tree8 plus the observation model in the E-step, the fission draw, the test maps and the
       scores.
     - Random-split controls: no call on either root (z_mid 0.0 / -1.5).
   - VERDICT (rnasep-misfit, ~02:00; `$H/handoff/rnasep-misfit/VERDICT.md`): RNase P's per-particle misfit is each particle's noise
     spectrum.
     - Residual power over 117-15 A varies with a true SD of 29 % (reliability .96). It is low-frequency and CTF-shaped, clusters
       by FoilHole and micrograph (ice), and follows from denoise_project's whole-band normalisation. Registration plays no part.
     - The fix, nvarB: a 6-band inverse-gamma prior fitted by marginal likelihood; the posterior mean enters the whitening and the
       fission draw. It is knob-free and collapses on IgM (prior CV <= .013).
     - Floor at r131: the balls outside the molecule go from 76 to 0.4 nats (z 0.3), and the rim controls from 245 to 29. A
       scalar per particle fixes 20 %.
     - Sparse pursuit with nvarB and the per-round control check: r131 best block .640 -> .70-.74 and span .694 -> .75-.77 (the
       same site in every run, ~(40, -33, 7) A). 692k 50k subsample: .737 -> .791 and .794 -> .853, beyond native on all 692k
       (.802 / .849). IgM arms within .003.
     - The per-round check is required: without the floor, the as-written stop runs away (held-out LL falls after 8 blocks).
     - Shared tree, 692k 50k subsample: the noise part of the root split is gone (.20 -> .07). The root now splits along
       amplitude (.38) into a 10 % low-amplitude, pose-unstable side (amplitude .52 vs 1.03; angular jump 20 vs 5 deg). Closed
       is not reached at level 1 (AUROC .495).
     - The user judged the low-amplitude group at the 692k root to be real RNA conformations, not junk (junk_screen=no). So that
       side is read as a candidate population: honest maps of both sides, and growth to depth 3.
     - r131 shared tree with honest y2 gate maps: running (job 55054551, root ~03:00).
   - r131 RESULT (~03:00): the shared tree is still one leaf under nvarB, so the hard gate is not met. The honest y2 map reads
     .81 / .62 on every combination (native or nvarB-relearned basis, base or nvarB E-step).
     - nvarB does remove the misfit mixture: the K = 2 XD held-out gain falls from .127 to .0095 nats per particle.
     - What remains are balanced proposals along the largest spread (47/53, thirds; d' 1.0-2.1), all no call.
     - Closed is a 7 % minority at d' ~.9 along one block's axis, too weak to form its own XD component at 44k.
     - So the port stays on hold. The next lever is the proposal: separability-ranked proposals (7.9 item 5), briefed to
       tree-latent (STATUS addendum 5) with r131, qtoy2's 1 | 2 pair and IgM a5 as targets, and every null as the guard.
   - 692k 50k subsample with nvarB, FINAL (misfit STATUS, "692k: the two sides of the root"; depth 3, 3.2 h on 36 cores):
     - The root split (K = 2, z 12.6, d' 3.98; amplitude .38, log-likelihood .29) is the only one taken. Both sides are leaves,
       and every proposal below reads closed at .48-.52.
     - Honest y2 maps, closed / open arm occupancy: all particles .31 / .85 (open); main side (45,121) .31 / .86; low-amplitude
       side (4,879, closed share .006) .21 / .51 (band-limited .19 / .20).
     - So the low-amplitude side holds less density in both arm regions, the open arm most. Its arm is weaker or elsewhere, not
       the closed arm. Described, not judged (junk_screen=no stands).
   - PORT STARTED (rnasep-misfit, worktree `SIMPLE_noisevar` on `5b17de846`): nvarB into flex's whitening, where flex already
     whitens per particle (`build%esig%sigma2_noise(shell, iptcl)`). It needs one residual pass at the node mean and 6
     two-parameter fits, ~30-40 s per 50k.
     - Default on, `SIMPLE_FLEX_NOISE_VAR=no` for A/B. refine3D is untouched.
     - Validation: native vs Python multipliers; native sparse_sites on r131; IgM and qtoyflex unchanged; the gate battery against
       pp_off with the r131 corner on both seeds. Committed only after the gates.
   - PORT DELIVERED, uncommitted.
     - Files: patch `$H/handoff/rnasep-misfit/noisevar_5b17de846.patch` (5 files, +345 lines, 2 comment lines); build
       `$H/wip_builds/noisevar_5b17de8`, fast gate 15/15; new module `simple_flex_noise_var` (in check_flex_dag).
     - Design: esig%sigma2_noise is per shell and half, not per particle, so the multiplier scales each particle's observation
       planes per band in `planes_batch_load`. Every flex consumer sees sigma2 x v; refine3D is untouched.
     - The master estimates v in prepare when vol1 is given; workers read `flex_pca_noise_var.txt`. It is on by default, and
       `SIMPLE_FLEX_NOISE_VAR=no` turns it off.
     - (a) Native vs Python multipliers: r .976-.993 per band. The rest is one per-particle factor, from flex's prep normalising
       each image by its background sd; with it removed, r .994-.999.
     - (b) Native sparse_sites at r131: controls 109.5 -> 9.4 nats, blocks 18 -> 9, span .711 -> .765, best block .659 -> .681.
       The native per-round check needs no change.
     - (c) IgM identical to the third decimal. qtoyflex closed .9005 -> .9010.
     - CAVEAT: on the toys the multipliers absorb a little state signal, because the residual is taken against the node mean
       only (qtoy2 root: v's holo AUROC .61 in band 1, CV .054). On RNase P, v carries no closed signal (.47-.53). The gates
       decide whether it matters; if they do, the residual is taken after the node's model, not the mean.
     - (d) Gate battery against pp_off ("nv" = the noise-variance build; "nvoff" = the same build with the switch off):
       - Nulls: one leaf each (qtoynull, olenull, micnull).
       - qtoyflex: passes (NMI .223 vs .225; closed .46 vs .47).
       - nvoff reproduces pp_off exactly on every fixture it ran.
       - qtoy2: NMI .756 vs .780, inside the seed spread (.707-.812). Closed and holo hold (.90 / .95). Combos 43 and 49 merge,
         and one node's second state becomes a pure junk leaf.
       - r131 seed 2 passes (corner .959 / .541 vs .960 / .544).
       - r131 seed 1 MISSES the corner gate by .003 (.947 / .543 vs pp_off 1.004 / .530). The miss is inside pp_off's own seed
         movement (.960 <-> 1.004), and nv's first r1 group recalls more closed rows (.60 vs .56; seed 2 .80 vs .62).
     - DECISION (coordinator): no commit as delivered.
       - Fixes requested: band edges in resolution, not box rings (rings 2-14 match 117-15 A only at box 36); a team-vs-serial
         tester for the two new parallel regions; default OFF (`SIMPLE_FLEX_NOISE_VAR=yes` enables); workers adopt only this
         run's file.
       - Then commit default-off after the fast gate.
       - Seeds 3-4 at r131 for both arms, plus nv qtoy2 seed 2. Default-on is the user's decision, with that 4-seed table.
     - COMMITTED default-off as `03506060d`: build 3 (`noisevar_e6b9ff844.patch`, md5 b030657d), with its own build and fast
       gate 15/15 (job 55066343, `$H/wip_builds/nv_03506060`).
       - It has the fixes: band edges scaled to the node's band limit (box 36 keeps today's rings), a team-vs-serial tester in
         `unit_heterogeneity` (suite `flex_noise_var`), a stale file retired so workers can't adopt it, and the switch
         `SIMPLE_FLEX_NOISE_VAR=yes`.
       - Builds 2 and 3 agree to 3e-7 on IgM.
       - SEED TABLE (r131 corner, best stage-B map, own data: closed >= .95 / open <= .6; closed-reference purity of the best
         depth-1 group is the independent readout, base share .069):

         | seed | pp_off map | nv map | nv gate | best-group purity (recall), pp_off / nv |
         | --- | --- | --- | --- | --- |
         | 1 | 1.004 / .530 | .947 / .543 | miss .003 | .078 (.56) / .083 (.60) |
         | 2 | .960 / .544 | .959 / .541 | pass | .085 (.62) / .106 (.80) |
         | 3 | 1.016 / .511 | 1.040 / .545 | pass | .081 (.54) / .116 (.80) |
         | 4 | 1.008 / .502 | .953 / .546 | pass (+.003) | .073 (.54) / .091 (.60) |

       - On the own-data map gate, pp_off passes 4 of 4 and nv 3 of 4. On the independent readout nv concentrates the closed
         reference rows more on every seed (mean purity .099 vs .079, recall .60-.80 vs .54-.62).
       - qtoy2 seed 2 with nv: NMI .778 (pp_off .780), closed .91 (.89), holo .94 (.94), and combos 43 / 49 separate again (.87 /
         .83). So seed 1's merge was that seed's draw.
       - The default-off build reproduces pp_off r131 seed 2 exactly.
       - Coordinator's recommendation: default ON. The closed state is not lost but better concentrated on every seed, and the
         one-seed own-data miss sits in the imprint-prone readout. The user decides.
     - Absorption check (nvarBG6, a signal-free background estimate on 6 bands): every 6-band or group model re-learns the
       closed block at .72-.74, so nvarB's closed gain is not signal absorbed into the noise level.
   - Checks under nvarB:
     - The fission copies' residual correlation across noise deciles goes from -.18..+.23 to -.02..+.04, and the fake-map null
       within noise terciles is flat (.02-.07).
     - Micrograph sigma2 groups remove 90-93 % of the outside floor (nvarB 99.5 %), and pooling over the 692k root adds
       nothing. A background-only estimate matches the residual at the same 2-band resolution, so the gain is the 6-band shape.
2. Data-derived band (sparse-ppca2). The IgM scan rows land, then the rule: gain above the solvent controls, per node or per
   block. Raw gain is out, since it follows the misfit floor on r131.
   - The rule is built into sparse_sites together with the profiled coordinates, validated on the six fixtures, and committed
     to wip after the fast gate.
   - That lets the site chain fire on IgM inside the tree.
   - IgM 10 % fix band scan (lp 15 lattice; round-1 best seed gain):
     - 25 A: 2,402 nats, all 5 arms in 4 sites;
     - 20 A: 1,081, all 5 arms, best readouts;
     - 17.5 A: 140, one arm weaker;
     - full band: 12.8, none.
     Rule (a), the cap with the largest leading gain, picks 25 A; every arm is recovered, with R^2 within .003 of 20 A.
   - The gain maximum sits at the scan's edge on every fixture: IgM at the coarsest cap, qtoy and r131 at the finest. So the
     range cannot place the peak yet. 19 runs on norm extend it: IgM fix and null at 30 / 35 A, and qtoyflex, qtoynull and
     r131 on a finer lattice reaching 12.5 A. The rule evaluation (a vs b) and the design note follow.
   - Profiled coordinates (`sparse_sites_profiled.patch`, build b4, fast gate 15/15).
     - The fit is unchanged (1e-4 of b3); 692k closed AUROC .802 -> .797 (best block), .849 -> .844 (span).
     - Every other readout is within .012, except the IgM 10 % null's moving arms, which lose up to .035 R^2.
     - Committed as `aac135dfe` on wip after the commit's own build and fast gate passed (15/15, job 55057869,
       `$H/wip_builds/prof_aac135d`).
   - BAND RULE (`$H/handoff/sparse-ppca/band_rule.md`, 39-run scan on diagnostic build b7d): rank caps by the PURSUED model's
     held-out gain, not the seeds'.
     - Seeds pick 31-36 A on IgM, where the six blocks fuse into one component and arm R^2 drops .01-.04. Detection wants coarse
       caps and the fitted model wants fine ones.
     - Pursued gain picks 17.7 A for the IgM fix (5 arm components, R^2 within .008 of the scan's best) and 20.4 A for the IgM
       null (no discreteness). qtoyflex and r131 get their full band, and qtoynull has 0 sites at every cap.
     - The rule: a quarter-octave ladder (2^(1/4)) from the node's band, down at most one octave, searched toward coarser caps,
       stopping once the pursued gain falls below the best so far after a site is kept. The constants are 2^(1/4), one octave
       and 2 SE.
     - Per-block bands wait for a per-block ring limit.
     - ACCEPTED, and briefed for the build (worktree `SIMPLE_band` on `5b17de846`): six-fixture validation, wall time per node,
       and r131 re-read once native nvarB lands.
3. Site chain (site-chain). v7 (sklearn-parity tied seed, b6, 15/15) worker tests run on norm. The b5 battery table goes into
   STATUS.
   - If v7 is better at 5 % with nulls at 0, it is the chain.
   - Once the band rule lands, the IgM chain trees rerun. If IgM gains in the tree and the gates pass, site_chain is committed
     to wip.
   - Then the factorial-output proposal (no build).
   - b5 battery done (`$H/handoff/site-chain/STATUS.md`, "Gate table", 14 rows). Every section-3 row passes and reproduces
     pp_on's partition:
     - qtoy2: 8 leaves, NMI .811, closed .89, holo .94;
     - r131: corner kept on both seeds;
     - nulls: 1 leaf each;
     - ole: 8 leaves.
     The chain never fired on a gate node, because the lp-15 sparse fit finds 0 blocks on every one. It stays unexercised until
     the band rule lands.
   - v7 is the chain (build b6, fast gate 15/15; worker tests 55057520; v7 columns in `round_table.md`).
     - Its seed now matches sklearn's tied fit (5 % a1 .180 against .181, a3 .162 against .164).
     - Why: sklearn's k-means++ start plus its 1e-3 stop keep the seed on the states. A fully converged tied fit leaves them
       (5 % a3 .164 -> .005).
     - 10 % fix: 5/5 arms DISCRETE, product .348 (v6 .321, lever 3 .338). a2 is accepted for the first time.
     - 5 % fix: 4/5 arms, product .148 (v6 .078, lever 3 .104). a2 (.119 at GAPw .29) is the only miss.
     - Nulls: 0/15 and 0/5.
     - No b6 battery: every gate firing had 0 blocks, so the b5 table stands. The battery reruns with the band rule.
   - Factorial output (`$H/handoff/site-chain/PROPOSAL_factorial_output.md`):
     - Per sparse site, it writes the accepted block with the lowest GAPw, once per lineage at the node that owns the sparse fit,
       in the existing factorial/ layout.
     - The tree's decisions are unchanged.
     - Approved as v8 on wip `a4886d747`, with no new key. It is judged by worker-test parity (.348 / .148) and empty tables on
       the nulls.
     - v8 (b9 on `5b17de846`, fast gate 15/15): parity is exact.
       - Table products .348 / .148, both null tables empty, and all 33 blocks' labels bitwise v7's, with the profile on or off
         (the worker never reaches the flex embedding).
       - Found on the way: `chain_blocks` read `nlines(coords.txt) - 1` rows, so the last particle was never seeded. It is a
         one-line fix (`site_chain_coords_read_fix.patch`), now built with v8 and rerun; that build is the commit candidate.
     - Commit candidate b10 (`site_chain_v8_plus_coordsfix.patch` on `5b17de846`): fast gate 15/15, products .349 / .145 with
       v7's accepted blocks, both null tables empty. It waits for the band rule and the IgM chain trees.
   - b11 = b10 + band3 (`$H/handoff/site-chain/STATUS.md`, 2026-11-08 afternoon):
     - Gate battery, 14 trees, all rows pass. Where the chain did not fire, the trees match pp_off exactly (qtoy2's 8 leaves at NMI
       .780; r131's maps to three decimals). r131 corner on own data: s02 .960 / .544, s01 1.004 / .530, and it holds at maxdepth 2.
     - The chain fired on 7 gate nodes and the band3 fit kept 0 sites each time, so the battery shows no harm but cannot show
       benefit. It never fires on r131, because every node there splits before it becomes a would-be leaf.
     - IgM stock twins: 10 % fix 2 leaves (combo NMI .034), 5 % fix 2 leaves (.001), nulls 1 leaf each. The 4 IgM chain trees are
       still running, and `sc11_igm_score.sbatch` scores chain against stock when they end. Cost per firing: sparse fit 2-15 min at
       32 threads, worker 2-5 s.
4. Pose-profile default (pose-tangent2). Seeds 3-5 per arm finish, then a recommendation: keep on, gate on the removed
   evidence share, or off. The committed default changes only with the user.
   - 5 seeds per arm (qtoyflex): the cost is systematic.
     - Closed F1 with the profile on is .315 +- .063, against .441 +- .012 off. Holo is .586 against .601, lower in every seed.
     - A K = 3 root appears on 0/5 seeds with it on, against 2/5 off.
     - Cause: the profile thins the root's site layer (fewer kept poles).
   - Gates: a 2-SE pose-share gate fails (qtoyflex R^2 .003 is about 14 SE). The truth-free evidence share separates cleanly
     (IgM .75-.87, elsewhere <= 0).
   - DECIDED (user, 2026-11-08: "Yes, default off"): `4afe93b80` turns it off by default, and `SIMPLE_FLEX_POSE_PROFILE=yes`
     turns it on. Gate trees compared with earlier pp_on arms must set the key.
5. Flex -> refinement priors (oracle, Part 7): the verdict when the chains finish.
   - VERDICT (igm-refs VERDICT.md Part 7; `results/part7/`): a fixed latent prior helps only where the seed already carries
     the arm.
     - 10 % fix: α = 1 beats lever 3 on all five arms at every round (+.012-.053), and the labels freeze after round 1.
     - 5 % fix: the old chain seed is a 1-3 % tail split, and α = 1 holds the chain there (-.011 to -.081 at round 3). It blocks
       the rounds that rescued these arms in Part 3.
     - Nulls: the prior does not stop round 1's carving.
     - The held-out α stays > 0 on every null block (.125-.25), so it is no null detector, and the site test stays GAPw. The
       held-out α is the safe weight: +.008-.049 over lever 3 at 10 %, within -.012 / +.004 at 5 %.
     - Caveat: β* sits at its grid top in 39 of 45 null and 5 % rounds.
   - Part 7b approved: the 10 % null and both 5 % fixtures, seeded from v7's sklearn tied refit, with no-prior / held-out α /
     α = 1 arms from that one seed, and a wider β grid.
   - Part 7b VERDICT (igm-refs VERDICT.md; `results/part7b/`): the fixed prior amplifies whatever split its seed carries.
     - From the tied refit seed it helps the 5 % arms: α = 1 minus no prior at round 3 is +.015 / +.021 / +.020 / +.003 on a1 /
       a3 / a4 / a2, and it holds near-true proportions where no prior pulls toward 50 %. Most of the gain over lever 3 comes
       from the seed itself.
     - It anchors and deepens balanced null carves. α = 1 deepens the 10 % null carve on all five blocks (+.010 to +.016 over no
       prior), so the null criterion fails.
     - The held-out α hits its grid floor in most rounds (SIMPLE α 51-68), where it freezes the seed. It is never the right
       weight.
     - RULE: the prior is only a polish, at α = 1, on a site GAPw has accepted, seeded by the tied refit. Never weight it by
       the held-out α, and never use it as an acceptance test.
     - The proposed GAPw check on the tied seed's null carves is already covered: site chain v7 accepts 0 / 15 and 0 / 5 on the
       nulls.
6. Every verdict goes into this guide and memory, pushed to fork wip. New jobs go to norm. Nothing goes to production.
7. The whole shared model (tree-latent; user: "what about the whole shared model").
   - Thread 1 decides RNase P. The model has two untested parts besides, and both run tonight.
   - Per-subtree refinement of the shared basis ("refinements per subset"):
     - At each accepted node, refit the pursuit on the node's members, on y2, since its stop decides structure.
     - Keep the parent's blocks (nested), and add blocks only by the knob-free stop. Then regrow the subtree.
     - Targets: qtoy2's conformer 1 | 2 pair (pp_on, which re-finds sites per node, gets NMI .811 against .766) and IgM a5's
       two 40 % states. Every null must gain no block.
     - If it resolves a miss with the nulls clean, it becomes the port's per-node step. Else the frozen root basis stands.
   - The knob-free band:
     - When sparse-ppca2's rule lands, the fixture battery reruns on the rule's bases. Today the IgM result rests on a fixed
       lp 20.
     - In the refined model a child's new blocks take the child's own band, so the frequency march lives inside one model.
   - Pose refinement per subset stays the port's job (stage A per product leaf).
   - Part 1 (`tl_nested.py`; at every node of the frozen K3-Sidak forests, parent blocks fixed, the root pursuit's settings and
     stop, on y2):
     - qtoyflex: no block on any of the 4 nodes (best seeds 0.1-1.3 nats, SE 0.1-2.7), so closed F1 .558 holds.
     - qtoynull root: no block. olenull and micnull have no root site.
     - qtoy2 r12 and r1: a seed passes round 1 (8.2 nats, SE 1.3; 5.0, SE 2.1), then decays under the group penalty. That is a
       power limit of y2.
     - So the pursuit also runs on y with the parents and nuisances profiled. Every y1 selection reads y1 only through the
       parent-and-nuisance span, so under white noise y's complement is untouched by selection. On real data this needs the
       per-particle noise model.
     - Leak checks: qtoyflex nodes, plus the IgM 10 % twin root and qtoynull root. The fitted-prior y variant is the contrast.
   - Part 2, FINAL (STATUS addendum 4; `math_reworks/tree_latent/refine/summary.txt`): no node gains a block, so the frozen
     root basis stands. PORT_DESIGN.md: "Per-subtree refinement of B: tested, not adopted".
     - Nodes tested: qtoyflex 4, qtoy2 15, IgM 10 % fix 11, the IgM twin root, and the three null roots.
     - Blocks added: 0 on y2, 0 on profiled y. The only fitted-prior y block (qtoy2 r1, seed 36, 20.4 nats) has no y2 or
       profiled counterpart, so it is the leak and is out.
     - Neither miss is basis-limited. qtoy2's seeds pass round 1, then decay under the group penalty; seed 36's gain lies in
       the parent span (30.6 fitted vs 2.1 profiled). IgM a5 stays below 2 SE. The limit is the proposal and the test.
     - At nodes the binding stop is the group penalty, not the 2 SE gate.
     - Trap: micnull's y2 root shows a 102-nat central seed (SE 5.0) that y and profiled y do not.
       - It is NOT the fission draw: the same y2 draw with the nuisances profiled has no such seed (best 5.6 nats, the seeds y
         shows).
       - CONFIRMED: micnull is synthetic and the draw already matched its spectrum (.985-1.006 per shell). A generator-spectrum
         draw keeps the seed (98.3 nats), while profiling the nuisances removes it.
       - The cause is the second-moment nuisance prior fitted on the noisier copy. Its centred low-frequency templates are
         under-absorbed (EM stops at its 200-iteration cap), and a central seed takes their variance.
       - RULE (PORT_DESIGN.md, fission step): anything fitted on a fission copy profiles the nuisances.
   - Separability-ranked proposals, part 1 (`tl_tree10.py`, `--pp 1`).
     - Design: a homoscedastic K = 2 XD per site on y1, i.e. a 1D unequal-weight split along the Fisher direction. Starts come
       from the 10/90 % tails along each principal and coordinate axis. Fitted on even members, ranked by odd held-out gain,
       proposed only above 2 SE. Sidak counts 2d tests per site.
     - qtoyflex: unchanged (closed F1 .558).
     - IgM 10 % twin: one leaf in all 5 trees. But on every arm the proposal finds a significant ~15 % minority (26-125 nats,
       d' 3-4.6), and 4 of 5 pass z_mid above the Sidak bar. sat5's side test (r at the continuum bound) rejects them all. So
       on continua the saturation test, not the multiplicity bar, prevents false splits.
     - IgM 10 % fix: product .321 (frozen .322). a5's 40/40 pair is still unresolved: a minority search does not split a
       balanced pair.
     - qtoy2 and r131 with nvarB still running. r131's closed is a 7 % minority, the case this proposal targets.
   - Separability-ranked proposals, FINAL (STATUS addendum 5; PORT_DESIGN.md: "tested, not adopted"). None of the three
     misses is resolved, and the must-holds pass.
     - qtoyflex and qtoy2: unchanged.
     - IgM 10 %: .321 (frozen .322).
     - r131 with nvarB: one leaf, and the honest gate fails (.81 / .62).
     - At each miss the best direction is right but hardens only a pure tail, and sat5 makes no call. qtoy2 r12: 877 particles,
       74 % conformer 2. IgM a5: 1,562, 97 % state 1. r131: 104 particles, closed .135 vs .069.
     - At r131 the most separable direction still tracks the residual misfit (log-likelihood .33). The K2 / K3 proposals no
       longer track misfit under nvarB, but carry no closed signal (.50-.51).
     - CONCLUSION: the states overlap along every single latent direction. A frozen-pose latent split cannot separate
       overlapping states as wholes, and a full-image round does.
       - Site chain v7: a5 .264 against the forest's .216, product .348 against .322, and a2 accepted.
       - The production tree keeps r131's corner through stage B and state_split.
     - DESIGN REVISION (done in PORT_DESIGN.md, section "The full-image round"): stage B is SEEDED by the latent proposals, not replaced by them. The forest's
       proposals and factorial states seed a full-image round on the node's subset, and an honest test (sat5 or GAPw) decides.
       This is the user's "refinements per subset are necessary".
       - The round takes the node subset (X1 members, stage-A poses), seeds as hard MAP labels on X1, and the seeding site's
         support as a mask. Part 7b (thread 5) supersedes the design note's soft-weight option: the prior only polishes an
         accepted site, at α = 1.
       - The test reads fission copies from the same draw under the per-particle noise model, or GAPw's own split. Seeds are
         tried in the forest's order, and Sidak counts the rounds tested.
       - It stays design-only until the site chain's r131 round reads honestly.
     - r131 ROUND TEST: NEGATIVE (site-chain VERDICT item 6; STATUS "RNase P r131"). The b10 worker ran at r131's root (44,334;
       closed share .069) on the native sparse fits, noise model on (9 blocks) and off (18 blocks).
       - Result: 0 DISCRETE blocks. Every block's groups are near-halves (17k-28k), and the purest group is 1.4-1.7x enriched
         (noise on: b7 g2, 17,149 rows, recall .64, purity .115).
       - Its maps read the node mean (.89-.90 / .59) against the .95 / .6 corner, so own data already fails.
       - Why: closed is a 7 % minority at d' ~.9 with no second mode, so the K = 2 seed halves the node and the round keeps that
         split. On IgM the arm states are 27-73 % groups, which the same seed and round resolve as wholes.
     - The production tree also halves r131 first ([22,484 | 21,850] at seed 2). It finds the corner one level down (r1 s02,
       .96-1.00), in a child with closed recall .80 / purity .106, the same enrichment as the chain's groups.
     - Next: the two-level test, on own data first. The chain runs inside its two best-enriched groups (sparse fit, seed, round,
       corner readout), and the logs show what pp_off's r1 stage B did (seeds, states, pose re-search or not).
       - If the corner appears, the design question for the user is an exploratory carve (not discrete) followed by an honest
         test one level down.
     - TWO-LEVEL TEST: NEGATIVE (site-chain STATUS, "The two-level path"). As child nodes, both enriched groups give sparse fits
       where no site survives (seed 29.8 / 33.9 nats keeps 3.6 / 3.0 once grown), so the chain has nothing to test.
     - WHERE THE PRODUCTION CORNER COMES FROM (`r131_ctrl/`, pp_off and nv seed 2): stage A's iterated full-image 2-state
       assignment, seeded by its fflex groups. The corner appears in stage A (.948 / .951), before stage B.
       - The tree's own r1_s2 labels, rebuilt at FIXED root poses, still read the corner (.972 / .963). Pose re-search and the
         focused PCG are not needed; the partition is the ingredient.
       - That partition is barely enriched in the reference closed rows (purity .085 / .106 against the share .069). The chain's
         purer b7 g2 (.115) reads .887. Reference-row purity does not predict the readout.
       - So iterated likelihood selection against closed-ish maps picks particles whose own data rebuilds the closed arm. The
         extremes verdict's fission re-read of the delivered tree's corner leaf gave .93 / .59 against 1.03 / .54 on own data.
       - KEY: the r131 corner gate, read on own data, largely measures selection imprint. Honestly read, no method shown so far
         passes .95 at r131: the shared model .81, the production corner leaf .93 (delivered tree), the chain's groups .89.
       - Tonight's own-data gate calls inherit this, e.g. the noise model's seed-1 .947 against the bar. This makes decision 9
         (the r131 gate's honest readout) the first decision for the user.
     - The r131 round test above was the planned decisive test (seeded from the noise-model sparse fit). It is negative because
       its K = 2 seed halves the node. The next test seeds the round from minority tails instead (thread 9).
   - Band rule build: qtoyflex picks its full band (9,480 nats, as in the scan) and qtoynull keeps 0 sites at all five caps.
     r131 hit an illegal instruction on a non-x6342 node and reruns pinned to x6342 (job 55059118). The IgM runs are going.
     - Cost: 1.5-8 min per toy node, 10-19 min per IgM node, 8-36 GB.
   - BAND LADDER DELIVERED (sparse-ppca2, 2026-11-08; `$H/handoff/sparse-ppca/band_ladder.md`):
     - Patch `sparse_sites_band3.patch` (md5 b5d9c457) on `7ae8ba4b1`; build `$H/wip_builds/band3_7ae8ba4`, fast gate 15/15; run
       jobs need `--constraint=x6342`. Not committed.
     - The v1 / v2 builds (`band_5b17de8`, `band_7ae8ba4`, `band2_7ae8ba4`) are superseded.
     - IgM 10 % fix: 18.2 A, 5 sites, one per arm (R^2 .523 .400 .485 .442 .400; the scan's best .523 .398 .485 .442 .401).
     - IgM 10 % null: 18.2 A, 5 sites, none discrete. The ladder has no 20.4 A step, and the gain is flat there.
     - IgM 5 % fix: 25.2 A, one site per arm (R^2 within .003 of the 20.4 A runs on four arms, a2 +.016).
     - qtoyflex and r131: full band, with the scan's readouts (closed .901 / holo .929; r131 best block .659, span .711, and
       .681 / .765 with the noise model).
     - qtoynull: 0 sites at every cap.
     - Cost: 1.6-2.3x one fit when a site is found, 4.0-4.7x when none is. 45-57 % of IgM time goes to bands that keep nothing.
       IgM nodes take 2.2-3.2 h at 39 GB. The 692k root would take ~3.5-4.3 h, or ~130 min with the first two bands as separate
       jobs.
     - Flags:
       - A dropped site is not undone, so the kept sites drift. This caused r131's 150-nat result at 18 A with the noise model
         on, and it predates the ladder.
       - Re-read r131's pick if the noise model becomes the default.
     - Next (thread 3, briefed 2026-11-08 morning): site chain b10 + band3 in one build (b11, job 55069844), then the IgM chain
       trees with stock twins and the gate battery.
     - User (2026-11-08): "not sure why its so much compute ... couldnt you bake this into the latent variable model itself ...
       we know how particles are related at this particular frequency and lets stick with that".
       - Why it costs so much: each rung is a fresh full pursuit, a node with nothing to find runs every rung, and about half of
         IgM's time goes to seeds that die while growing.
       - Briefed to sparse-ppca2: a band learned inside one fit, either a per-site per-shell prior variance by marginal
         likelihood (ARD over shells) or the per-shell half-set agreement of the modes as Wiener weights. It is compared with
         the ladder on the same seven runs, at a target cost of ~1 fit.
       - The ladder stays uncommitted until that reports.
       - IN-MODEL RESULT (sparse-ppca2, 2026-11-08 afternoon; `$H/handoff/sparse-ppca/inband_band.md`): none of four versions
         replaces the ladder. Each loses the IgM 10 % fix entirely (no site), where the ladder finds all five arms.
         - Versions:
           - v1: a per-site per-shell prior.
           - v2 / v3: a per-block band edge by held-out gain. v3 is batched: r131 .678 / .743 (.740 / .772 with the noise model),
             qtoyflex .900 / .929, at 1.0-1.5x.
           - v4: one node-wide per-shell weight from the even/odd agreement of the node covariance. On IgM it falls below one
             half at 110-330 A, far coarser than the arms' 18-30 A.
         - Mechanism: a seed needs its band fixed before its first iteration.
           - Learned during the fit, the band never cuts in time. v1's edge stayed at ring 43 while the seed's energy fell 14x,
             and v3's narrowing collapsed seeds faster.
           - A node-wide probe finds a global low-frequency direction, not the local arm motions.
         - Next (briefed 2026-11-08): the one design left open. Measure the per-shell even/odd covariance agreement inside the
           seed balls and pool it into one band before phase 1.
           - If it loses the IgM fix, the ladder is final. Add the two cost changes (caps 1 and 2 as concurrent jobs; a round-1
             null screen at box 36), confirm the picks, build, run the fast gate and commit to wip.
         - Flag unchanged in every version: the dropped-round drift (v3 r131 with the noise model: LL/ptcl .244 -> .170 over 12
           rounds).
         - v5 (the seed-ball band) LOST on IgM. Its pooled agreement inside the seed balls stays .50-.76 on every IgM shell
           (likely each particle's contrast times the mean map), so it never cuts the band. Round-1 seeds gained 8.6 / 4.7 nats,
           as at the full band. qtoyflex fell to .844 / .875.
         - THE LADDER IS FINAL. Commit `f4be24169` on branch `sparse-sites-band-ladder` (off wip `20e972dfd`, worktree
           `/Users/atifao/SIMPLE_ladder`): band3's rule, plus two cost changes.
           - Cap 2 runs as its own queued sub-job while the node runs cap 1. It writes to its own subdirectory with its ring and
             random state; the parent picks by the unchanged rule, and runs cap 2 in-process if the sub-job fails.
           - A round-1 screen tests caps 3+ with the node band's seeds, only when that band's round 1 finds no seed. IgM never
             triggers it.
           - Build and fast gate: job 55086887. Then the seven ladder fixtures must reproduce every pick, block count and gain.
             After that it goes on wip.
           - Cap jobs need x6342 nodes (USE_ARCHOPT builds). The parent sets `SBATCH_CONSTRAINT=x6342`, and tree launches need
             the same. SIMPLE's template mail line (`--mail-type=FAIL --mail-user=nil`) reaches nobody.
           - Revised to `eb1d99ede`: cap 2's side job had crashed on restoring the saved random state, so earlier checks
             passed only through the in-process fallback. It passes the fast gate (15/15).
             - All four box-36 fixtures reproduce the ladder exactly with cap 2 as its own job. qtoyflex's wall drops from 297
               to 176 s.
             - The IgM and r131 checks (55087717-19) are running.
           - `eb1d99ede` must NOT go on wip. On the IgM 10 % check cap 2 won, but BeeGFS refuses the folder rename, so its
             sites were never taken over and both runs reported 0 sites. The fix is in `e4fbef54e` (on wip `cc4af8a10`):
             build 55088875, then six validation runs.
           - On `e4fbef54e` the box-36 fixtures reproduce exactly. The IgM 10 % fix now takes over cap 2's sites and reproduces
             the ladder (ring 36, 8 blocks in 5 sites, R^2 .523 / .400 / .485 / .442 / .400) in 1 h 23 instead of 2 h 13.
             A log-count fix makes the final commit `6aa8492cb`; it is building, with reruns after its gate.
           - ON WIP as `2f41cbc76` (rebased `6aa8492cb`; fast gate 15/15). All seven fixtures match band3's picks and readouts.
             The IgM 10 % fix runs in 1 h 23 instead of 2 h 13, the null in 2 h 03 instead of 3 h 10.
             - Untested: ladder plus the new noise-on default on IgM / qtoy (the validation ran noise off). It is opt-in
               through sparse_sites=yes, so wip's defaults are unchanged.
8. Honest r131 gate (decision 9; extremes-repro, briefed 2026-11-08):
   - One fixed fission pair of the r131 project (X1 / X2), with the noise drawn at sigma2 x v (the per-particle band-wise noise
     model), checked with the fake-map null.
   - Each gate tree runs on X1. The corner leaf is the one whose X1 map reads most closed, and it is read on X2.
   - The bar is the 95th percentile of matched-N random tails and random subsets from the corner's parent node, read on X2.
   - First use: pp_off vs the noise model on seeds 1-2 (3-4 if close), which decides the noise-variance default. Then one script
     for the battery's r131 arm.
   - SEEDS 1-2 RESULT (`$H/handoff/extremes/GATE.md`; battery script `scripts/ext_gate_tree.sbatch`). Corner, own X1 -> honest
     X2, against the bar:
     - pp_off seed 1: r1_s1 (10,459), 1.074 / .599 -> .893 / .641 against .892: PASS (+.001).
     - noise model seed 1: r2_g1_2 (21,354; r2 unsplit), .882 / .534 -> .843 / .607 against .856: FAIL (-.013).
     - pp_off seed 2: r1_s1 (10,983), 1.004 / .569 -> .870 / .607 against .883: FAIL (-.013).
     - noise model seed 2: r2_s2 (12,620), .919 / .563 -> .864 / .620 against .855: PASS (+.009).
     - Every corner loses .04-.18 of closed occupancy when rebuilt from X2. Read honestly, the arms split one-all and each pass
       is narrow.
     - The stricter null is the latent tail along the parent's own X1 latent, which reads as closed as the corner about 5 % of
       the time.
     - Root splits reproduce within an arm (79 % / 77 % agreement across seeds) but not across arms (50-56 %). The noise model
       moves the root split onto a different axis.
     - FINAL, 4 seeds (`$H/handoff/extremes/GATE.md`):
       - Seeds 3-4: pp_off .873 against .877 and .862 against .871 (both fail); noise model .857 against .861 (fail) and
         .851 against .844 (pass, but the corner sits at the base closed share .071).
       - Honest passes: pp_off 1 / 4, noise model 2 / 4. Own-data passes on the old gate: 4 / 4 and 0 / 4.
       - Not separated: every call is within about one null sd of its bar, and the corners' closed excess over their parent
         is the same in both arms (+.005).
       - The noise model shrinks the own-data imprint: its corners lose .04-.06 from X1 to X2, pp_off's lose .11-.18.
       - The r131 corner gate cannot decide the noise default; pp_off's own-data lead was imprint.
       - Battery script `scripts/ext_gate_tree.sbatch` (ARM, BUILD, SEED; ARMENV or ARMKEYS exported). Cost 209.7 core-h
         plus 10.1 for the pair.
9. Making the shared model work (user, 2026-11-08: "subsetting is necessary (with non random seeds derived from the latent) and
   refinement is necessary. I only know this empirically"). Both tests approved 2026-11-08 (user: "Yes, run both").
   - Why, from mixture-learning theory:
     - A 2-class fit at r131 prefers to halve the largest continuous spread. A 7 % minority is therefore a local optimum,
       reached only from a start inside its basin: the latent seeds.
     - Provable mixture algorithms initialise in a low-dimensional projection and then iterate in the full space (Kumar &
       Kannan 2010). The latent is the projection. Refinement (free 3D references, all bands, poses) is the full-space
       iteration that recovers what the projection lost.
     - Subsetting raises the minority's share of the variance, pi(1-pi) (7 % -> 30 % triples it), and removes competing spread.
     - Full-space iterations on the same data also create imprint, hence the honest test.
   - Overnight, each RNase P failure lacked one ingredient:
     - the shared tree had no full-space refinement;
     - the chain's r131 round started from balanced K = 2 halves;
     - the production corner has all three (the root halving with closed recall .80, fflex seeds, iterated full-image
       assignment) but is read on its own data.
   - Design:
     - The latent is the proposal layer. The shared basis is kept, with minority tail proposals beside the XD K = 2/3.
     - Each proposal seeds a full-image round on the node subset: hard MAP labels on X1, the site's support mask, fixed poses
       then local search, and the noise model in the likelihood.
     - An honest fission test against matched random groups accepts; sat5 / GAPw type the result (discrete state or
       reproducible extreme).
     - Exploratory carves act as lookahead and are merged back unless a child passes honestly.
   - Decisive test (proposed): at r131, a full-image round on X1 from the best latent tails, from random halves and from the
     production fflex groups, read on X2 against matched random tails, with the IgM nulls as the false-call guard.
   - APPROVED and briefed (2026-11-08):
     - Carve test (rnasep-misfit; `$H/handoff/rnasep-misfit/carves/`):
       - The shared model with the noise model on all 692k carves along each node's leading proposal wherever no honest split
         passes, to depth 4.
       - Readouts: reference closed share / recall per node, and honest y2 reads against matched random subsets, compared with
         the production lineage (thread 10).
       - Pass: a lineage reaches ~7 % closed share at recall >= ~.3 by its own carves, with its honest map above the bar.
       - Guards: the IgM 10 % null twin, qtoynull, qtoyflex.
     - Tail-seed refinement test (extremes-repro; `$H/handoff/extremes/`):
       - At r131 on the honest gate's X1 / X2 pair, production stage A on X1 from three kinds of seed: latent tails (picked with
         the closed reference, so a mechanism test), random seeds, and the fflex groups.
       - Each is read on X2 before and after the round, against the random bar. A fixed-pose variant runs for the best tail and
         one random seed.
       - Guards: the IgM 10 % null twin and fix.
     - CARVE TEST SO FAR (2026-11-08 afternoon; `CARVES.md`):
       - Guards:
         - qtoyflex holds: closed best-leaf F1 .571 against .558. Only the root split honestly; 11 nodes were carved.
         - qtoynull holds trivially: no site survives its basis.
         - The IgM 10 % null twin FAILS once. S3:r11 at level 2 split honestly (K = 3 pair 1,176 / 2,617, z_mid 3.3 against the
           Sidak bar 3.32) under a carved "not discrete" parent. That is 1 pass in 73 tested nodes (177 tests); all five
           random-split controls read no call.
       - The per-node Sidak bar does not control the tree once carving multiplies the nodes tested.
         - Briefed: re-read every pass against a size-and-depth Bonferroni split, alpha_node = alpha (n_node / n_root) /
           (maxdepth + 1), then Sidak within the node. The shares come from the X1 tree.
         - That bounds the false-pass rate over the whole tree.
         - RESULT (a re-read of the recorded tests, `bar_reread.py`):
           - Under the tree-wide bar the IgM null's pass goes away: z_mid 3.33 against 4.11 at node share .216. The null
             tested 74 nodes, not 73.
           - qtoyflex's root still passes (3.98 against 3.75), so its F1 .571 stands.
           - The r692 root and all 6 honest splits through level 2 pass both bars. The tightest is r23: z_mid 5.93 against
             4.62.
       - r692 (job 55076637): the root split honestly, K = 3 into 601k / 82k / 8.8k. The tree is on level 3, and the table with
         both bars comes at about 18:30-19:30 EST.
       - CARVE TREE RESULT (from `carves/tree_r692/tree.log`; tree done 18:01; rnasep-misfit's formal table to follow). FAILS
         the pass bar (~7 % closed at recall >= .3).
         - Main lineage: r 1.26 % -> r1 (638k) 1.31 -> r12 (642k) 1.31 -> r121 (596k) 1.27 -> leaf r1212 (248k) 1.94 %,
           recall .55. That only matches production's depth-1 node r1 (252k, 1.94 %), three levels later.
         - Best leaf: r1223 (14k), 3.46 % at recall .055.
         - Production for comparison: r131 6.94 % at recall .35 by depth 3; r1312 8.3 %; r2322 9.4 %.
         - The honest root split (z 42) does not concentrate closed. Then levels go to peeling off small groups (r11 8.6k, r122
           50k) before the bulk is carved, and maxlev 4 stops at production's depth-1 equivalent.
         - Honest maps: r1223 reads .44 against a random maximum of .35, and r1212 .34 against .29. The closed signal is real
           but thin.
         - Formal table in `$H/handoff/rnasep-misfit/carves/`. The tree-wide bar refuses only r221's honest split. The job's
           remaining step adds a lookahead for each production node.
       - USER (2026-11-08): "yeah dont think we can get out of refitting at each node". The shared-basis tree is closed as a
         replacement for the tree's per-node refits on RNase P. Production's edge is its K-way root and per-node latent refits,
         which make closed lead where the tree cuts. It is not refinement (pose test) and not one global basis (this tree).
     - REFINEMENT TEST SO FAR (r131, honest pair):
       - Every round ends in a near-half split (about 19k of 44k), whatever the seed. The round's state stays near the base
         closed share (.069).
         - Random seeds (rtail1, rtail2, rsub1, rsub2) clear the bar 3 of 4 times, at closed share .062-.077.
         - Latent tails (pc2lo, site2hi, blockhi) clear it 1 of 3 times, at .076-.089.
         - Site-masked blockhi clears it 0 of 1 times, at .098.
         - rsub2 has both halves at the base share (.069 / .070) and still reads .857 against .853. Clearing the bar after a
           round does not need any closed signal in the seed.
       - IgM guard: a 97.5 %-pure 5 % a1 tail grew into a 44 % state, 64.5 % pure (a1 NMI .068 -> .051), still far above its
         bar.
       - The seed decides the axis. Rounds from different seeds land on mostly unrelated splits (50-66 % agreement). pc2lo,
         rtail2 and site2hi share one axis (64-74 %), near the noise-model trees' roots.
       - The fflex round reproduces the tree's own root split (94.3 % agreement with pp_off seed 2). It is the most
         closed-enriched round: 23,381 particles, share .085, recall .64, honest .862 against .852. No latent tail beat it.
       - Coordinator's read: the round as built (production stage A, K = 2) pulls any seed to about half the node, so a
         minority seed does not survive it.
         - CONFIRMED in code: neither stage-A path has a class-size prior.
           - The production round (prob_neigh) labels by a bare argmin of the distance over states
             (`eul_prob_tab%assign_rows`).
           - The fixed-pose round (state_eval) labels by an argmax, with a log-prior only if `SIMPLE_STATE_PRIOR` is set.
         - Trajectories: every coherent seed keeps 80-97 % of its particles while its state grows about x4 in iteration 1, to
           ~40 % by iteration 3 and 43-48 % by iteration 6, on r131 and IgM alike.
           - The fixed-pose round absorbs as much, so the pull comes from the assignment rule, not pose refinement.
           - The trajectories cannot separate the missing prior (~2 nats per particle at 7 %) from the small state's noisier
             reference.
         - Briefed: one fixed-pose round with `SIMPLE_STATE_PRIOR` = log class fractions re-estimated each iteration, from the
           IgM a1 tail and from r131's site2hi.
         - PRIOR RESULT, r131 site2hi (seed 2,217, closed .091):
           - With the EM prior the state holds at 1,936-3,105 by iteration 6, at closed .083-.092. The zero-prior control on
             the same harness and score grows to 20,427 (.081), exactly as the earlier round did.
           - So the missing prior drives the absorption.
           - With the prior, only 39 % of the seed's members remain: the group keeps its size while its members churn. It
             reads .910 honestly, against a bar of .927.
           - Stage A needs a class prior, but at r131 the prior alone does not concentrate closed.
         - IgM a1 with the prior ends at 29,371 (purity .651, NMI .068 -> .031). a1's true state is half the node, so this
           waits on its zero-prior control.
         - FULL DESIGN ROUND, EM prior plus the seed's own site mask at fixed poses (`$H/handoff/extremes/ROUNDS.md`; 294.4
           core-h for all rounds):
           - IgM a1: 30,619 at .875 purity, keeping 96 % of the seed, a1 NMI .068 -> .210. Production's round gives .051,
             the prior alone .031, and the site chain's own EM round .448. It was still growing at the last iteration.
           - r131 site2hi with site 2's mask: it holds 1.8-2.4k but churns at the seed's closed share (.091), honest .918
             against .931.
           - r131 blockhi with block 2's mask: 5,046 at closed share .128 (2x base, recall .21, 70 % of the seed kept). Its
             maps read below random sets of the same size (.883 against .891; bar .926).
           - My first brief paired site2hi's seed with block 2's mask. The two share no voxel (87 A apart), and that round
             collapsed. A seed must use its own site's mask.
           - IgM guards: the unmasked round fails both; the masked fix round passes (.095).
         - VERDICT: at r131 no full-image round, from any latent seed, concentrates closed in a way the honest maps can see.
           - The class prior is needed, and the mask makes the round work on IgM, but neither finds RNase P's closed state.
           - Open question: at r131 the closed reference rows and the arm-occupancy maps disagree (blockhi as seed and after
             its round). That needs explaining before closed share can judge rounds there.
           - The carve test at 692k (r692 table) is the remaining half of this thread.
10. Why the current tree finds closed where it does (user, 2026-11-08: "I think we have to figure out why we get the close states
    at the nodes we do in the current tree"; rnasep-misfit, briefed 2026-11-08; notes `$H/handoff/rnasep-misfit/lineage/`).
    - Tree: the 692k integration tree (recipe17d, `$H/integration/f692_r17d_s1`).
    - Closed paths: root -> r1 -> r13 -> r131 -> r1311 -> r1311_s2 / r1311_s1, r131 -> r1312, and r2 -> ... -> r2322. Own-data
      closed readouts along the main path: .45 -> .67 -> .87 -> .90 -> 1.00.
    - Per split: the route and seed source; what each side tracked (reference closed share / recall, and the AUROC of noise level
      v, amplitude, defocus, micrograph / FoilHole, view, pose change and the parent's latent dimensions); and the closed AUROC of
      each node's own latent by depth.
    - The puzzle: r131 holds ~7 % reference closed rows (~1.3 % at the root), but the corner r1311_s2 holds 816 of 12,010 (6.8 %),
      no more than r131, while its own-data map reads 1.00. Honest reads of r131's leaves come from thread 8's fission pair.
    - RESULT (rnasep-misfit, `$H/handoff/rnasep-misfit/lineage/LINEAGE.md`; memory r17d-closed-lineage-attribution-2026-11-08):
      - Closed share by depth: r 1.26 % -> r1 1.94 -> r13 4.11 -> r131 6.94 -> r1311 5.57 (r1312 8.34) -> r1311_s2 6.79. On the
        other lineage: r2 1.31 -> r23 2.92 -> r232 5.71 -> r2322 9.37.
      - The concentrating steps are binary splits along ONE dimension of the node's own flex latent: r13 -> r131 (z1 AUROC .82)
        and r23 -> r232 (z1), each taking 93 % of the parent's closed rows; then r232 -> r2322.
      - The K-way levels above select good particles (SIMPLE corr / amplitude .59-.60, poses that moved less), with some closed
        enrichment.
      - Noise level, defocus, micrograph and FoilHole are flat at every split (.47-.53); views / 2D classes read .53-.58.
      - Refitting the latent per subset adds .02-.04 closed AUROC and never beats the root latent's .786. The closed direction is
        already in the root latent, and the tree carves along it where it leads.
      - The corner:
        - r131 -> r1311 went against the closed rows.
        - r1311 -> s2 cuts along a dimension the closed rows barely follow (closed AUROC .56; 4.2 -> 6.8 %). It is the weakest
          split that passed the re-test (z 4.87 vs 3.71), and stage B moved 16,009 of 22,464 particles away from the latent's
          groups.
        - So its 1.00 fits the imprint pattern; the honest reads decide.
        - By the reference rows, the most closed leaves are r2322 (9.4 %) and r1312 (8.3 %).
      - Coordinator's read for the shared model (thread 9):
        - The tree concentrates closed by carving along a closed-carrying latent dimension at the node where that dimension
          leads, whether or not the cut is discrete, after K-way levels that keep good particles.
        - The shared model's honest discreteness test refuses such carves, and its 692k proposals (50k subsample) never led with
          closed (.48-.52). So the carves, and the conditioning that makes closed lead, look like the missing pieces rather than
          a better basis.
        - Whether refinement adds closed is not shown here; thread 9's test answers it.
      - Running (approved): one map per lineage node at entry poses, then noise-model sparse sites on r13, r1311, r1312, r23,
        r232 and r2322 (~20 core-h). It asks whether a shared-basis site carries the closed direction where the tree cuts.
      - SITES RESULT (`$H/handoff/rnasep-misfit/lineage/LINEAGE_SITES.md`): native sparse sites with the noise model on rebuilt
        node means, 4-16 min per node.
        - r13 (13 blocks): the leading site coordinate ranks the tree's cut to r131 at .859 (the tree's z1 .821, latent LDA .847).
          Closed: block 1 .754, span .847 (z1 .718).
        - r23 (12 blocks): the cut to r232 sits on block 2 (.753 / .778; span .834 against z1 .819), not on the leading
          coordinate (.501). Block 2 reads closed at .836 against z1's .743.
        - r131 (9 blocks): closed block 1 .681, span .765. The cut to r1311, which went against the closed rows, is on block 2
          (.740).
        - r1311 (where the corner cut is made), r1312, r232 and r2322 (8-22k): 0 blocks. Every grown block adds 0-29 nats held
          out.
        - Read: the shared basis carries both carves that bring closed in, from its own coordinates and more sharply than the
          tree's z1. It has nothing at the deep nodes, which matches their weak re-tests (z 4.87 / 5.33) and the corner's imprint
          pattern.
        - For the carve rule: at r23 the carve is not on the leading coordinate.
11. Full 692k tree on wip defaults vs production (user, 2026-11-08: "Yes, run it tonight"). Session wip692; notes
    `$H/handoff/wip692/STATUS.md`; run `$H/integration/f692_wip_s1`.
    - CANCELLED before launch (user, ~17:00: "I dont want a safety check ... I dont want to run it if its not going to meaningfully
      change results"). The launcher stays prepared in wip692. No full tree runs until a lever shows a clear gain on node replays.
    - Correction to the baseline below: lf has SEVEN closed leaves at closed-arm >= .6 (r1312, r1321, r5111, r2112, r5121,
      r1322, r1311, all `_g1_2`; the 2026-11-03 union run used all seven). Three is only the count the user judged by eye on
      Nov 2.
    - Build: the noise-default flip on wip `0a29cfbee` (rnasep-misfit; fast gate and smoke check first), with the OpenBLAS
      variant matched to `build_lf`.
    - Keys and launcher exactly as `f692_lf.sbatch` (recipe17d keys, seed 1, all 692,392 rows). Expected 620-870 core-h, 5-6 h.
    - One run covers everything since 2026-11-02 (master merge, depth_wiener, noise model). It is compared with production's last
      692k run `f692_lf_s1` and with `f692_r17d_s1`.
    - Scorecard: production's own (leaf scorecard plus docking landscape; N = wip, O = lf). Rules fixed BEFORE the run:
      - Gate (hard): at least one leaf the user judges closed. Losing the closed state rejects the defaults, whatever else
        improves.
      - Same: every row inside the span of lf and r17d, the only yardstick of run-to-run variation at 692k. Purest closed leaf
        .095 / .094; best closed-leaf F1 .115 / .119.
      - Better: the closed state kept AND a closed row clearly past both, with the landscape no worse (all nine 2 A swing bins
        filled, best arm CC median >= .79) and wall within 7 h. "Clearly past" means purest closed leaf >= .11, best closed-leaf
        F1 >= .14, or 3-leaf union F1 >= .16 (lf .137).
      - Worse: fewer user-accepted closed leaves, purest closed leaf < .08, best closed-leaf F1 < .10, a shrinking landscape,
        or wall past 7 h.
    - It is one seed, and the maps are the tree's own data, so any call is provisional. The closed reference rows lead. A clear
      result needs a second seed to confirm it, and a 692k twin on the production tip (`merge_master2`) to attribute it.
12. Master's multi-state development (user, 2026-11-08: "look on master for new developments regarding multi state ... I think
    poses matter a lot but if multi state is scrambling membership too much it may be detrimental").
    - origin/master has 9 commits since our merge base `6634307d9`. One concerns states: `cac010e40` (2026-11-08), "multistate
      solve3D: draw the state from the posterior, final assignment at beta 1". Report:
      `doc/implementation_notes/completed/multistate_posterior_state_draw_report.md` on master.
      - The state is one draw from the tempered joint posterior over state and direction, p ~ exp(-T (L - L_min)), in place of
        the closest state. Labels stay hard.
      - Multi-state runs end with a final assignment at prob_beta = 1. The fill's budget and the annealing stall are per state.
        nstates = 1 is bit for bit unchanged.
      - There is NO population (class-size) term: every state has the same prior, as in our stage A.
    - Its own benchmark (20k simulated ribosome assemblies, 30 / 60 / 10 %) bears out the user's worry:
      - The 10 % state is never isolated, and the 60 % state splits over two.
      - At the low-resolution stages (22.6 / 20.2 A) about half the labels change every iteration, and the state posterior
        spreads over 2.2-2.7 of the 3 states.
      - Labels only organise at the end of stage 4 and in the 6 A stage.
    - For the tree:
      - Do not merge `cac010e40` into stage A/B without the gates. Its draw re-assigns membership even more than an argmin,
        it has no class prior, and a tree-wide prob_beta = 1 is the one arm on record that lost the r131 corner. The final
        assignment at beta 1 is not that arm, but it needs the same check.
      - It agrees with thread 9: multi-state re-assignment at low resolution scrambles membership and drains minorities.
    - Lever that follows (proposed): refine poses per group with membership HELD (single-state refine3D on each group). Compare
      at r131, on the honest pair, against the same groups at the parent's poses and against the multi-state round from the same
      seed.
      - Groups: the latent's label-free tail, blockhi (upper bound), the fflex group, and matched random groups.
      - Readout: the validated contrast readout.
      - Poses refined on X1, maps on X2.
    - APPROVED 2026-11-08 (user: "Yes, queue it"), briefed to extremes-repro.
      - Groups: blockhi seed (2,217), the blockhi prior+mask minority (5,046), a label-free ~5k tail of the top block's
        leading coordinate, and the fflex closed-side group (23,381).
      - Arms: (a) the parent's poses; (b) single-state refine3D on X1 from the parent's poses (stage B keys and iterations),
        maps on X2. The multi-state round from the same seed is the third arm.
      - The random bar gets the same single-state refinement (8-12 groups per N), because refinement alone may sharpen any
        subset. IgM a1 tail plus 4 random groups is the guard.
      - Readout: c = closed - open, applied only after the re-score's positive control shows the readout sees the reference
        closed rows.
      - Arm (c), added at the user's request ("freeze the memberships and just refine the poses, or only allow a very small
        number of state reassignments ... injecting the flex weights into the refinement"):
        - After (b), a fixed-pose state pass (state_eval) with SIMPLE_STATE_PRIOR = log class fraction plus an anchor: log
          w_flex for the fflex group, a one-hot on the current state for the others.
        - A particle moves only at posterior >= .95.
        - Two alternations of pose refinement and damped re-assignment.
        - Unlike the 11-06 flex->refine test (calibrated alpha, no class prior, inside freely re-assigning rounds; it switched
          itself off at r131), membership here starts frozen and can only move on strong evidence.
    - RE-SCORE, FIRST PASS (no new maps; contrast c = closed - open against the same random bars):
      - Gate corners: c passes one corner per arm (pp_off seed 3, noise model seed 2), with three calls flipped. The arms
        stay unseparated.
      - Seeds: only pc2lo clears its bar. blockhi (.184 closed) reads +.315 against .369.
      - Rounds: fflex, pc2lo, rtail1 and rtail2 keep their passes; rsub2 (random) loses its pass.
      - Within comparable groups the reference share and c rank alike: pose-refined rounds rho .79, fixed-pose rounds .99, gate
        corners .90; the seeds at N 2,217 only .16. In the rounds the open arm carries much of the signal (rho -.59 / -.93).
      - So labels and maps agree on ranking. But enrichment to .08-.13 closed does not move a map past its random bar.
      - Running: the positive control (r131's 3,075 reference closed rows on X2; optimistic, since the rows were picked on X)
        and a rigid-core term for blur.
    - POSITIVE CONTROL PASSED: the reference closed group reads c +.487 against a bar of +.359. The two-arm readout sees a pure
      closed group.
    - POSE TEST VERDICT (`$H/handoff/extremes/POSETEST.md`, `tables/posetest.md`; 76 core-h): at r131, membership limits the
      closed signal, not poses.
      - Arm (b), membership frozen and poses refined alone: every 2-5k group loses contrast like random groups of its size.
        blockhi -.062 against random -.081 +- .022; bpm1 -.069 against -.082 +- .008; blk1top -.079 against -.082; combined
        p = .15. At 23k (fflex) nothing moves.
        - The shift does not grow with closed share (28 groups: Spearman -.17, slope +.04 +- .10).
        - Even IgM's 97 %-pure a1 tail loses arm amplitude when refined alone (172 -> 154). Solo refinement of a 5k subset
          degrades its map.
      - Arm (c), damped re-assignment (EM class prior + anchor, move at posterior >= .95, two passes): closed shares stay
        flat (blockhi .184 -> .184, bpm1 .128 -> .131, blk1top .091 -> .093). Every change in c falls inside the random
        groups' changes.
        - Design flaw: for a 5 % group the class prior (-2.94 nats) cancels the .95 margin, so members leave at even odds
          while a newcomer needs 5.9 nats. Small groups shrink: random 21-30 %, blockhi 33 %, groups along real variability
          4-8 %.
        - Flex weights as the anchor are weak (median .50 nats, 90th percentile 1.87, under the 2.94 margin). Only 125 and
          then 30 of 44,334 particles move, so they just freeze the partition (closed .085, recall .645).
      - IgM guard passes: the a1 tail grows 5,000 -> 5,354 at .97 purity; random groups shrink to ~4,700 at .50.
      - Read: at r131, refining poses per subset or re-assigning with a damped rule adds no closed signal beyond what the same
        step gives a random group. The groups' membership (at most .18 closed, ~.09 label-free) is the limit. Raising it
        needs more per-particle evidence or a better membership rule, not refinement.
13. Per-node levers on production's closed lineage (user, 2026-11-08: "Yes, run the replays"; session flex-refine-flex-weights,
    notes `$H/handoff/node-levers/STATUS.md`).
    - Premise, from the user: per-node refits cannot be avoided, so the levers act inside each node's refit.
    - 2 x 2 factorial:
      - the noise model in flex, on (wip default) or off (SIMPLE_FLEX_NOISE_VAR=no);
      - a class prior in the state assignment, on or off. This is the log of each state's current fraction, added to the score
        in eul_prob_tab's assign_rows and the neighbourhood table's state choice. It sits behind an env switch, default off,
        and equal fractions reproduce the argmin.
    - Nodes: r17d's r131, r13, r1, then the root. Each is replayed once from production's saved node through the tree's own
      calls (f2r harness, maxdepth 0). The off/off arm on the new build is the baseline.
    - Readouts: the node latent's closed AUROC; each child's closed share, recall and concentration factor against
      production's own children (1.94 / 4.11 / 6.94 / 8.34 %); at r131 also the honest contrast on the fission pair.
    - Priced before the root (ask if over ~600 core-h). The class prior goes to the gate battery only if it wins.
    - Class-prior build `84d4bb75d`: fast gate 15/15. Production's pruning had removed the parents' recB0 halves and recB maps,
      so they are regenerated for r13, r1 and r (recB maps rebuilt from the stage-B project, a documented stand-in).
    - INCIDENT 19:13: the second replay harness reused a root holding a symlink into production. It wrote 1,599 self-links
      into f692_r17d_s1, and its rm -rf deleted the contents of production's node_r13/recB0 and node_r13/recB.
      - Restored: recB0/b.simple from a 19:10 copy, and both hetero_step.sh. The coordinator verified 0 links left.
      - Lost: the step logs, exit files, node_rows bookkeeping, and the pcg_support and RESOLUTION files in those two
        directories. No result, split, map, site or coordinate file was touched.
      - The harness now refuses symlinked or reused destinations and never removes through a link.
      - Production is NOT read-only: the session's permissions refused the recursive chmod. Making the tree read-only is left
        to the user (`chmod -R a-w .../f692_r17d_s1/1_hetero_tree`, reversible with `u+w`).
      - The first four r131 arms ran with the parent's mask; they are kept as a variant, and the corrected arms run in fresh
        roots.
    - r131 RESULT (corrected arms, production's eight-state mask):
      | arm | noise | prior | best latent dim | children closed share | closed child | split |
      | --- | --- | --- | --- | --- | --- | --- |
      | production | r17d | no | .699 | .056 / .083 | x1.20, recall .209 | z 13.1 |
      | a | off | off | .698 | .087 / .051 | x1.25, recall .226 | z 13.8 |
      | b | on | off | .713 | .069 / .070 | x1.00, recall .187 | z 14.9 |
      | c | off | on | .698 | one leaf (12 particles left in state 2) | none | fail |
      | d | on | on | .713 | one leaf (4 left) | none | trivial |
      - The baseline reproduces production.
      - The noise model sharpens the latent (.698 -> .713), but the node then splits on an axis unrelated to closed (x1.00).
        The noise-model gate trees showed the same root-axis shift.
      - The class prior as specified (MAP with current fractions) kills the 2.6 % seed: -3.6 nats per particle from
        iteration 1, so one leaf. It is not the lever.
      - FLAG: the noise model is default-on on wip since this afternoon. If r13 / r1 also lose closed alignment, revert it
        to opt-in.
    - r13 RESULT (closed share .0411):
      - production .718, x1.69 / recall .353, z 17.2;
      - baseline .722, x1.74 / .355, z 27.4;
      - noise model .732, x1.71 / .354, z 27.9;
      - class prior: one leaf (29 particles left; 22 with both levers).
      - The noise model does NOT repeat r131's loss here. In r131's parent-mask variant the arms were also equal (x1.36 /
        x1.36). So the noise model shows no stable effect on closed alignment, either way, and the default stays pending r1.
      - The MAP class prior removes the split at both nodes: the seed minority is charged log(pi_min / pi_max) (3.6 nats at
        r131) from the first iteration. A floored or tempered prior is not run: a 1/K floor is no prior for K = 2, and the
        EM-prior rounds already showed a held minority with flat closed share.
      - Cost ~5 core-h per r131 arm, ~9 per r13 arm. r1's four arms are running. The root stays on hold unless r1 shows a gain.
    - r1 RESULT (closed share .0194; 23-30 core-h per arm):
      - Production: 6 sites, site 2 dominant (margin 16.0 vs 4-5.6), 3 children, closed child r13 x2.12 / recall .380.
      - Baseline a (off/off): 4 sites at a near-tie (9.482 vs 9.466), site 1 seeds, 2-way split 132k / 120k with NO
        concentration (x1.10, recall .295). The new build's baseline does NOT reproduce production at r1.
      - Noise model b: 5 sites (9.862 vs 9.675), site 4 seeds, 3-way split, closed child x1.90 / recall .389, i.e.
        production's structure. The gain comes through the seed-site choice flipping at a near-tie.
      - Class prior: a 4,867-particle state, not enriched, or one leaf.
      - HARD-GATE CONCERN: the wip build loses production's r1 -> r13 concentration in its baseline. Possible causes: the
        build's site layer, or r1's parent maps (regenerated stand-ins for the root).
      - Briefed:
        - arm a on production's own builds (r17d's or build_lf, and merge_master2) with the same inputs;
        - arms a / b at a second SIMPLE_SEED.
        About 100 core-h; the root stays on hold.
      - CHECKS (55091355-8, ~25-30 core-h each): NOT a wip regression.
        - Build check, arm a on the same inputs:
          - production's own r17d binary (sha256 = the manifest's): 2-way split, x1.08, two sites (8.95 / 5.57);
          - merge_master2 = the wip baseline exactly: x1.10, sites 9.48 / 9.47 / 7.34 / 5.45.
          - The latent is identical across the three builds (.705). Every build misses production's dominant site-2 pole on
            these inputs; the only stand-in is the parent's recB map.
        - Seed check (new build, SIMPLE_SEED=2): baseline finds the 3-way split, x2.00 (recall .395); noise model x1.88.
          - With seed 1 (baseline x1.10, noise x1.90): the baseline's miss was a coin flip at a near-tie.
          - The noise model hit both seeds without beating the baseline's hit (x2.00) or production (x2.12). No gain; at most
            a robustness hint (2 / 2 against 1 / 2).
        - Reading: r1's closed concentration rests on the site layer's pole choice at near-tie margins. None of six replays
          found a pole as dominant as production's (16.0 vs 4-5.6). This is a reproducibility weakness of node site
          discovery and ranking, the same on r17d and on wip.
    - VERDICT (r131, r13, r1; ~490 core-h): neither lever concentrates closed at any node. The MAP class prior removes the
      split, and the noise model gives no gain. The root is not run. The open issue is reliability: the r1 closed split is a
      near-tie seeding outcome.
14. Robustness of the tree's decisions (user, 2026-11-08: "If this is not replicable, the methods must be redesigned ...
    considering the model itself how can we make it more robust").
    - Reliability test at r1, briefed: arms a / b at SIMPLE_SEED 3-6, 8 runs, ~225 core-h. HIT = a child at >= 1.5x r1's
      closed share. So far a 1 / 2, b 2 / 2.
    - Diagnosis: every split rests on hard, greedy, single-choice decisions taken from one noisy statistic:
      - the seed site (r1: 9.48 vs 9.47);
      - one site's poles as the seeds;
      - K;
      - accept / reject.
      Particles never move between branches afterwards. An error at a near-tie near the top removes the closed lineage, and
      nothing later can repair it.
    - Redesign options, cheapest first:
      1. Uncertainty-aware seed choice. Score each site on half-sets or seeds and compare the top site's margin with its own
         spread. Seed as now only when the top site is clearly ahead; otherwise use option 2. This is stability selection
         applied to the split choice; it adds no knob beyond the replicate count.
      2. Joint seeding over tied sites (the user's 7.5 idea). A mixture in the joint coordinate of the tied sites' focused
         flexes, with K by held-out likelihood and the map-level honest test deciding. It removes the coin flip
         structurally: the closed direction can't be dropped because it ranked second.
      3. Lookahead at uncertain nodes. Develop the top two candidate splits one level and keep the one whose children carry
         the stronger reproducible (honest) structure. Costs x2, at uncertain nodes only.
      4. Not now: a global re-assignment over all leaves at the end. Today's results (re-assignment dilutes, the MAP prior
         drains) argue against it until 1-3 are in.
    - Test: the same r1 seeds as the reliability test. Target: hit at every seed with nulls at one leaf. Then the closed subtree
      (r1 to leaves) at three seeds, reproducible closed leaves being the publishable claim.
    - User (2026-11-08): "if there are cheap tests for this go ahead. but remember we are trying to produce the most diverse set
      of conformations possible". A fix must ADD states at a near-tie, not pick one.
    - Cheap test briefed (f2r session, Python on the replay outputs, < 30 core-h, no Fortran):
      1. Are r1's tied sites different motions? Their masks, overlap, and pole-map differences across all replays.
      2. A joint-coordinate mixture over the tied sites at r1's entry poses, K <= 6 by cross-half held-out likelihood.
         Readouts: closed share and recall per group (1.5x rule), and the number of honestly distinct group maps (arm
         occupancy, swing).
      - The question: does the joint mixture hold the closed group at the seed where single-site seeding missed, and give
        more distinct conformations?
Waiting for the user:
- decision 1 (depth_wiener default);
- decision 3 (state_eval port; not needed by the chain);
- decision 6 (fission split test: no as it stands);
- decision 9 (the r131 gate's honest readout): DECIDED 2026-11-08, honest re-read (section 8 item 9; thread 8);
- the held precision fix;
- a production candidate;
- the noise-variance default (flex): committed off as `03506060d`. Seed table in 7.10 thread 1: own-data corner 3 / 4 vs 4 / 4, but closed-reference purity higher on 4 / 4 (.099 vs .079). Recommended: on. Honest re-read (thread 8, 4 seeds): not separated (honest passes 1 / 4 off vs 2 / 4 on, all within ~1 null sd); the noise model shrinks the own-data imprint (.04-.06 vs .11-.18). DECIDED 2026-11-08 (user: "Yes, default on"): on by default on wip, `SIMPLE_FLEX_NOISE_VAR=no` turns it off. ON WIP as `938ed9a8d` (cherry-picked from `dd02eff12`). The build `$H/wip_builds/nvon_dd02eff12` passed the fast gate 15/15 and links build_lf's OpenBLAS pair. The smoke runs on qtoynull, with the variable unset, show the master estimating the noise variance in both flex steps and all 4 workers adopting it. The ladder's r131 pick under the noise model is already re-read: ring 15, 9 blocks, span .765.

- tree-latent (shared model, RNase P tests):
  - Runs on the r131 and 692k sparse fits, plus the port design note.
  - If it keeps the r131 corner on honest maps, or isolates the closed lineage at 692k without tracking nuisances, the Fortran
    port starts (user: "port it if it's useful"). Otherwise the failure is written up.
- site-chain:
  - Variants of the round: no round, and refine3D's state_eval round. The 11 b3 trees run, then the section-3 battery.
  - A variant that passes the IgM worker tests and the gates is committed to `hetero-tree-wip` only, never to `hetero-tree`.
    If it needs the state_eval port (decision 3), that is flagged in the morning summary.
- site-coupling (factorial vs joint) and the 692k sparse PPCA result are folded in when they land.
- M-estimation retry (user: "maybe try M-estimation since these responsibilities seem more rigorous now"), queued with
  tree-latent after the RNase P tests.
  - Why it failed before: M-est added nothing because the tree's posteriors were one-hot or null-level, and soft EM at fixed
    poses noise-locked through full-image E-steps.
  - What changes now: tree-latent's responsibilities come from the deconvolved site-space mixture (posterior precision A,
    memberships on y1).
  - The test: each leaf's map built on y2 from hard members vs from responsibility weights. It is one M-step, so the E-step
    stays in latent space. Optionally a Student-t robust weight from the latent Mahalanobis distance.
  - Scored: against truth maps per state, minorities first, and the r131 closed-arm readout on hard vs soft maps.
- site-chain: the first chain-exercising null, `sc3md1_micnull`, is one leaf (sparse PPCA found no site there).
- M-est interim (tree-latent, qtoyflex, oracle references): soft weights trade bias for noise.
  - Closed leaf: deviation FSC .512 -> .611, error .863 -> .668, map FSC .736 -> .822, but the closed-specific slope falls
    .509 -> .235.
  - Holo leaf: correlation and FSC up, error and slope down. Arm leaf: slightly better on all four. Robust is close to soft.
  - Why: weights from y1 with maps from y2 give a confusion-weighted mixture of the state means.
  - Asked for: an UNMIXED estimator, mu = C^-T M with the mixture's confusion matrix, for soft and hard maps alike.
- tree-latent on RNase P (STATUS addendum 2; PORT_DESIGN.md, forest form, ON HOLD): FAILS the hard gate on both roots.
  - r131: no call, 1 leaf; honest y2 maps read closed .81 / open .62.
  - 692k: K = 3 [577k / 99k / 16k] with both pairs DISCRETE (z 31 / 35), but closed shares .014 / .008 / .006 (AUROC .48-.50).
    The splits track per-particle amplitude (.21-.46) and log-likelihood (.30-.36): misfit, not the closed state.
  - Why: fission cannot remove a per-particle property, the nuisances were integrated under a fitted prior, and a
    two-population amplitude difference saturates the group maps.
  - Cost at 692k: E-step 24 min (53.8 GB); root proposals 23 min on a 50k subsample; one sat5 test 51 s.
  - One misfit-removed rerun approved (amplitude + nuisances profiled; sat5 groups within amplitude x log-likelihood strata).
    If the best 692k split still tracks misfit, the real-data misfit model becomes the prerequisite, not the tree.
- DECISIVE RERUN (STATUS addendum 3; VERDICT.md 2026-11-08 section): amplitude + nuisances profiled, sat5 groups in 3x3
  amplitude x log-likelihood strata.
  - r131: 1 leaf; the honest gate fails (closed .81 / open .62).
  - 692k: K2 [616.5k / 75.9k] DISCRETE (z 20.3, was 43.2), but closed AUROC .483. The split tracks log-likelihood .42 and
    amplitude .30.
  - Coordinates + amplitude separate closed at LDA .745 / .856, so the signal exists; the splits just do not find it.
  - Mechanism: a per-particle noise level off the model correlates y1 and y2 by ~e/(2+e), so misfit survives fission.
  - STOP: the real-data noise model (rnasep-misfit) is the prerequisite. It must enter the whitening AND the fission draw.
- K = 3 final on qtoy2: closed .910 (pp_on .891), holo .962 (.944), NMI .766; the 1 | 2 pair is still merged.
- M-est final:
  - Hard maps are the default. Soft dilutes every state leaf; robust is the compromise and beats hard on large IgM leaves.
  - Unmixing restores the slope but over-corrects by up to 38 % and blows up on small leaves.
- M-est unmixing restores the slope at a noise cost: qtoyflex closed .51 -> .91 (error .86 -> 1.98); qtoy2 closed .84 -> 1.00
  (error .19 -> .24). IgM 10 % follows.
- Flex -> refinement priors revisited (user: "perhaps we should also revisit flex --> refinement posteriors"; oracle session,
  igm-refs Part 7).
  - Why it may differ from the 11-06 verdict: on IgM stage B never resolves the arms, so a prior has nothing to fade against;
    and the posteriors are now sparse-site ones with precision A, not the miscalibrated whole-molecule latent.
  - The test: state_eval with SIMPLE_STATE_PRIOR = alpha log w from the tied seeds' responsibilities, against lever 3's forced
    seed, with alpha by checkerboard held-out and at 1. Rounds 1-3; null carving, and alpha on the nulls.
- Idea (user: "sparse ppca kind of reads like a nonuniform filter"): both choose a local model per region by noise-calibrated,
  held-out criteria on a spatial grid. The NU filter picks a local band for the MEAN (candidate bank + Potts smoothing); sparse
  PPCA picks local supports for the VARIABILITY (a group lasso, calibrated by sign flips, + the per-shell prior). A unified
  "nonuniform sparse PPCA" would let each site's modes take their own local band from `simple_nu_filter`'s machinery.
  Untested; for the paper's framing.
- Nothing is promoted to production (`hetero-tree`) without the user. Every verdict goes into this guide and memory, pushed to
  the fork.
- BAND (user: "really the band should be data derived"):
  - Sparse PPCA finds all IgM arms at lp 20 and none at the tree's LPF 15. Giving lp 15 the lp-20 seed size changes nothing,
    so seed geometry is ruled out.
  - The band cap must come from the data: the cap that maximises the leading seeds' held-out gain per node, or per block in
    the nonuniform-filter style. sparse-ppca2 scans caps from 25 to 12.5 A on IgM, qtoy and r131 and proposes the rule.
  - CAUSE CONFIRMED: lp 15 with only the band cut at 20.4 A recovers all 5 IgM arms (matching lp 20 within .007; round-1
    gain 1,081 vs 12.8 nats). The 15-20 A shells carry no arm signal but most of each seed's cells, so the shared shrinkage
    prior pulls every seed to zero.
  - Box-36 scan: qtoyflex is best at the full band, and qtoynull finds no site at any band. r131's gain is largest at the full
    band, but its closed state reads best at 21-26 A, and its solvent floor grows ~3x with the band. So raw gain alone
    (rule a) follows the misfit floor on real data; the rule must measure gain above the solvent controls. IgM scan pending.
- USER DIRECTION (late evening): "We need to pursue this to work on the RNase P. If it works on IgM like regimes that is a great
  result." So the shared tree model stays the main methods thread. The RNase P blocker is per-particle misfit, worked on two ways:
  - tree-latent: the misfit-removed rerun (amplitude and nuisances profiled, sat5 groups within amplitude x log-likelihood
    strata), on r131 and 692k;
  - rnasep-misfit (new session; `$H/handoff/rnasep-misfit/`): what drives the misfit (residual band powers, amplitude and
    log-likelihood against defocus, micrograph, the 2D-class registration's rotation / sub-pixel shift, class, part), and
    observation-model fixes, each scored on three things:
    - the r131 solvent-control floor (109.5 nats now, which should fall toward 0);
    - whether the latent stops tracking amplitude and log-likelihood;
    - the sparse blocks' closed AUROC.
    The fixes: a per-particle band scale, a per-particle noise variance, and an interpolation-aware noise spectrum.
  - rnasep-misfit step 1 (r131 root; `$M/rnasep_misfit/tabs/`): the per-particle misfit is the NOISE LEVEL, not signal.
    - Band-wise residual power relative to the model's white per-shell noise has a true SD of .29 over 117-15 A, with
      split-half reliability .96. Nuisances and sites remove ~2 % of it, and it carries no closed signal (AUROC .48-.51).
    - It is low-frequency and anti-correlated with 6-3 A power (r -.62): denoise_project normalised every image by its
      whole-band background sd.
    - It clusters by FoilHole and micrograph (LOO r .42-.70), rises with particles per micrograph (neighbours in the box) and with
      defocus. Near the first CTF zero it follows the particle's own CTF signal (r +.72), so noise shaped like the transfer looks
      like variance to every local mode. That is the ~100-nat floor.
    - The 2D registration plays no part (|r| <= .08), so the interpolation-aware fix is dropped.
    - The tree-latent "log-likelihood" tracks |site coords| and the in-band noise.
    - Amplitude tracks SIMPLE's contrast (r .89) and CARRIES closed signal (AUROC .655), so profiling amplitude away costs closed
      information. tree-latent has been told.
    - Lead fix: a per-particle band-wise noise variance by empirical Bayes (inverse-gamma prior per band, posterior mean in the
      whitening). It is scored on the solvent floor (12 rim + 12 corner controls).
    - Added checks: the noise estimated outside the dilated molecular mask vs from the residual, so it cannot absorb closed
      signal; closed AUROC under every fix; and per-micrograph sigma2 groups as the cheap intermediate. SIMPLE's sigma2 groups on
      these stacks are the 10 stack parts.
    - Port entry if it wins: a per-particle (or per-micrograph) sigma2 per shell in flex's observation operator.

## 8. Decisions waiting for the user

1. depth_wiener on by default: gated and cleared, but no truth-map benefit (section 6), so the measured case is for `depth_wiener=no`.
2. Drop the chart-coordinate port (bundle ref `wip/tree-ports-charts`) unless the seam-free coordinate passes 7.2's rule (it does not, on
   r131 or at 692k). The CTF-flip trap (section 9) does not change this (|CTF| rerun: <= .003).
3. Keep state_eval as a developer harness: add it to `hetero-tree-wip` (after a build and fast gate) or leave it in the bundle.
4. site_split: after gate 1 and qtoyflex.
5. Of the three ideas, the flex prior, depth wiener and the global coordinate all failed their tangible-benefit tests on the
   hard fixture (section 6). What remains is the site layer itself: arm-shaped masks from the data (T1b's `self` masks are running),
   focused coordinates beside the node's, and latent-built references assigned per site. Port only after the IgM fixture and its
   null twin, the 10 % fixture and the section-3 gates.
6. The split test: replace `state_split`'s locality guard and effect floor by a fission test (labels from X1, maps from X2,
   tau 0.6; 7.8 fission VERDICT: the null passes, qtoy2 keeps its states, .02 accuracy cost). Needs a real-data noise model
   (per-group spectrum; r131's stacks are non-stationary) before the gates; the band split is the fallback. Prototype (user:
   "Yes, gate nulls first"; 7.8 `fission-split`):
   - With the guards off, fission merges qtoynull and micnull (self-selection) and keeps olenull's contrast split, which the
     locality guard merges.
   - The effect floor is redundant.
   - The candidate rule is fission + the locality guard (VERDICT in 7.8 `fission-split`). It passes every null and keeps
     qtoy2's states (NMI .717 vs .811), but r131 fails the hard gate at tau 0.6 AND at tau 1: its root, a near-tie unrelated to
     closed, never splits honestly, so the depth-1 corner is never reached.
   - As it stands, fission cannot replace the guards without losing RNase P's closed corner. Today's guards stay.
   - Options for the user: keep today's guards; recover power by cross-fitting both copies (data thinning); or accept that the
     honest test makes r131's root a no-call and judge the closed lineage at 692k, where it has more particles.
   Fission removes
   noise-made splits, not carves of a real continuum (the IgM null twins' root split follows the arm angles); igm-refs Part 5's
   saturation test is the continuum check that would complete the rule.
7. Pose profile in flex's embedding: `5f80720e0` on `hetero-tree-wip` after the gates passed (7.8).
   - DECIDED 2026-11-08 (user: "Yes, default off"): off by default since `4afe93b80`; `SIMPLE_FLEX_POSE_PROFILE=yes` turns it on.
   - Why: over five qtoyflex seeds it costs closed (F1 .315 against .441), and no tree-level gain has been shown.
   - The flex-precision fix (stacked patch, 7.8) is correct in principle and closer to truth at the root latent, but neutral to
     slightly negative on the gates (qtoy2 closed merged on 2 of 3 seeds, 0 of 3 without it): HELD. The deconvolution's input
     is wrong today, but the gate trees do better with it as it is.
8. Site layer: sparse PPCA's pursuit (knob-free stop: 2 SE + solvent controls, leave-one-out, support components) + one masked
   round, the best data-only chain (7.8 igm-refs Part 3).
   - Acceptance: GAPw on every block (Parts 4 / 4b: all IgM arms, no null or extra block).
   - The fission saturation test serves as a check for curved continua, not a replacement: with round-1 labels it loses 3-state
     arms, and it fixes qtoyflex only with true groups.
   - Neither test isolates qtoyflex's closed state from realistic labels. That is a seeding problem, which tree-latent's XD
     splits and the tree's own splits must cover.
   - INTEGRATION (user: "Yes, integrate it", "keep integrating and testing autonomously"; site-chain session, worktree
     `/Users/atifao/SIMPLE_site_chain` at `6206971c8`, notes to `$H/handoff/site-chain/`): an opt-in site layer running sparse
     sites -> per-block K = 2 seed from the block coordinates -> one masked round -> GAPw per block. It is add-only, like
     `site_split`. It runs on the IgM fixtures first, then the gate battery.
     Design (approved 21:00):
     - Key `site_chain=yes` (default no; exclusive with `site_split`). It replaces only the would-be-leaf worker's site source
       (sparse PPCA's blocks instead of `discover_sites` / `discover_islands`) behind a new step S_SITESPARSE.
     - Per block: the seed is a K = 2 GMM on the block's coordinates, with several starts; one round of `site_split`'s
       classifier (cross-half references, niter 1); then GAPw. The lowest-GAPw DISCRETE block gives 2 children.
     - The sparse fit is learned once per lineage and reused down the subtree (the shared-basis principle). lp = LPF 15.
     - Extra chain-exercising arms: the three nulls at maxdepth 1 (the false-accept test), r131 at maxdepth 2 on both seeds
       (hard gate), qtoy2 at maxdepth 2.
     - The factorial split (every accepted block at once) waits for site-coupling's verdict.
     Build (21:30): `$M/site_chain/b1` = `6206971c8` + `site_chain_v2.patch` (6 files, ~300 lines; ports recipe on x6342 with
     OpenBLAS, matching the pp_on / b3 baselines), fast gate 15/15 (job 55039984). The approved changes are in: best-of-5 seeds,
     and a would-be leaf reuses its own or its nearest ancestor's sparse fit.
     IgM trees running (site_chain=yes, igm keys, maxdepth 4): `$M/gates/sc1_igm{10fix,10null,fix,null}` (55040098-101). The
     gate battery follows. The stock same-run tree comes from collapsing the chain splits (`sc_collapse.py`).
     - BUG found on the way (also hits `site_split=yes`): when stage B leaves a single state, the node called decide() and was
       then forced to S_DONE. That dropped a submitted per-site split and left its rows unplaced. COMMITTED as `b712021b2` (one
       line; the step stays on S_SITESPLIT, as the S_SPLIT branch already does). Audit: the six earlier `site_split=yes` trees made
       9 per-site submissions, all at 2-state nodes and all resolved, so no earlier run was hit.
     Build b2 (`site_chain_full_v3.patch` = the fix + `site_chain_on_sdone_fix.patch`): fast gate 15/15 (55040431).
     - Children skip blocks already used on their path; seeds take the best of 8 starts; lp 15.
     - 13 trees run on b2 (21:24, `$M/gates`):
       - IgM: `sc2_{igm10fix,igm10null,igmfix,igmnull}` (55040553-6), plus stock same-binary twins `sc2off_{igm10fix,igm10null}`.
       - Chain-exercising arms: `sc2md1_{qtoynull,olenull,micnull}`, `sc2md2_{qtoy2,r131,r131_s1}` and `sc2_qtoyflex`
         (55040559-65).
     - The section-3 battery follows the first clean firing.
     - b2's first firing crashed: the worker's command line began `dir=` instead of `prg=` (the key-order trap). Build b3
       (`site_chain_full_v4.patch`) fixes it, passes the fast gate 15/15 (55041807), and runs the 11 chain trees
       (55041962-72).
     Worker-only validation (55041863; the b3 worker on the native lp-20 root fits, against truth), MIXED. Lever 3's Python
     round is in brackets.
     - 10 % fix: a1 / a3 / a4 accepted, with NMI at or above lever 3 (.448 [.377], .356 [.349], .329 [.318]). The 3-state arms
       a2 / a5 are lost (GAPw .33-.39 [.19-.22]).
     - 10 % null: 1 / 15 false accept (b9, GAPw .199) [0 / 13].
     - 5 % fix: a1 lost (no call); a3 / a4 accepted on almost nothing (NMI .05-.06).
     - 5 % null: 0 / 5.
     - Cause: site_split's one cross-half round re-randomises weak seeds (47-50 % of labels change at 5 %, 10-28 % at 10 %),
       where lever 3's state_eval round strengthened them (5 % a3 .012 -> .139). This is the "labels fall short" gap of 7.1,
       now isolated from discovery and masks.
     - Next, in parallel: the same tests with the round skipped, and with refine3D's own state_eval round from the bundle ref
       `wip/state-eval`, in a separate patch. Porting state_eval into wip is the user's decision 3.
     The seed is the bottleneck, more than the round.
     - Seed-only (no round): 10 % seeds a1 .389, a3 .145, a4 .138; at 5 % every seed is <= .013.
     - sklearn on the same coordinates: the full-covariance GMM picks unbalanced high-likelihood optima at 5 %, while a TIED
       covariance cannot. Tied lifts all five 5 % arms to .10-.18 and matches full at 10 %.
     - Full seeds, Fortran round vs state_eval round (port builds, 15/15; `state_eval_port.patch`): the Fortran round is better
       at 10 % (a3 .356 vs .276), state_eval lifts weak 5 % seeds about twice as far, and neither fixes a dead seed. state_eval's
       maps also remove the 10 % null false accept (b9 .97).
     - Tied seed + Fortran round (v6, build b5, fast gate 15/15; worker tests 55045233) is the winner:
       - 10 % fix: 4/5 arms DISCRETE (a1 .448, a3 .414, a4 .382, a5 .267; a2 not discrete), product .321 (lever 3 .338).
       - 5 % fix: 3/5 (a1 .154, a4 .087, a5 .067), product .078 (lever 3 .104).
       - Nulls 0/15 and 0/5; b9 is gone (1.21).
     - Tied seed + state_eval round: running (55046338).
     Trees so far:
     - Chain-exercising nulls (maxdepth 1) are 1 leaf each.
     - qtoy2 maxdepth 2: 9 leaves, NMI .811 = pp_on, closed .89 / holo .94.
     - On all 7 gate-node firings the lp-15 sparse fit found 0 blocks, so the gates only show the chain does no harm.
     - r131 maxdepth 2 never fires (both depth-1 nodes split), so it gives the hard-gate readout only.
     - Stock twins (b2): the 10 % null is 1 leaf; the 10 % fix root splits 52.5k / 47.5k into two leaves, so the chain fires at
       r1 / r2.
     - On b5: IgM chain trees `sc5_*` (55046502-5) and the section-3 battery `sc5_*` (55046506-12).
     ROUND TABLE (`$H/handoff/site-chain/round_table.md`; 33 blocks x 6 variants; the lp-20 root fits): the tied seed carries the
     gain, and the round matters little once the seed is tied. Accepted arm blocks (product of their labels) / null false
     accepts:

     | seed + round | 10 % fix | 5 % fix | 10 % null | 5 % null |
     | --- | --- | --- | --- | --- |
     | full, no round | 3 (.160) | 0 | 0 | 0 |
     | full + Fortran | 3 (.264) | 2 (.031) | 1 (b9) | 0 |
     | full + state_eval | 3 (.215) | 2 (.057) | 0 | 0 |
     | tied + Fortran | 4 (.321) | 3 (.078) | 0 | 0 |
     | tied + state_eval | 4 (.280) | 3 (.081) | 0 | 0 |

     - The b9 false accept exists only with the full seed + Fortran round.
     - a2 is never DISCRETE. The tied seed loses 5 % a3, which only full seeds catch.
     - The chain as delivered is tied + Fortran (v6, b5). DECISION 3 IS NOT NEEDED by the chain; the state_eval port
       (`$M/site_chain/state_eval_port.patch`, 5 files, b3se 15/15) stays a separate patch for the user.
     - Remaining lever: the Fortran tied fit converges worse than sklearn's at 5 % (seed a1 .061 vs .181, a3 .005 vs .164).
       Asked for convergence parity (k-means++, sklearn's n_init / tolerance) as v7.
     tree-latent profiling cost: profiling amplitude out of the site coordinates moves their closed AUROC only .706 -> .700
     (r131) and .849 -> .844 (692k). Amplitude's own closed signal (.655 / .703) never entered the coordinates.
     Nuisance check (site-coupling trap):
     - The seed coordinates (sparse PPCA's coords.txt) INTEGRATE the 12 nuisance latents under their fitted Gaussian prior
       (`simple_flex_sites_fit.f90:426-449`); the port's design note said "profiled". The round profiles (amplitude only, flat
       prior).
     - Within one site's 2-D seed this is not the cross-site trap, but it may blur seeds, and it matters for any factorial
       output.
     - Asked of the sparse-ppca2 session: profiled coordinates beside the integrated ones, validated on the fixtures and r131 /
       692k, as a follow-up patch to `b3ced0021`.
     SIDE FINDING to check: qtoyflex closed F1 .159 (holo .595) on this base's stock tree, against .47 for ports4 / ports7 at the
     same maxdepth, with no chain split involved. The pose profile (5f80720e0, default on) was never gated on qtoyflex. An A/B
     (profile on vs off, one binary, 2 seeds) has been asked of the pose-tangent2 session.
     A/B RESULT (`$H/handoff/pose-tangent/VERDICT.md`, qtoyflex section; build_poseprof, no chain code):
     - Closed F1, on vs off: seed 1 .159 vs .474, seed 2 .413 vs .433. Holo ~equal. So the .159 is the profile at seed 1, not
       the chain.
     - Mechanism, a NEAR-TIE FLIP: the profile is a ~2 % perturbation of the root flex states. At seed 1 the root's pole rule
       then keeps one site, so the root splits K = 2 instead of K = 3, and r1 (81 % of closed + 24 % of holo) is merged by the
       P/N 2 effect floor at 1.45. Off kept r1 at 2.66.
     - The profile does not remove closed signal: the closed map is 0.6 % in the tangent span, and the root latent's closed
       AUROC is .746 vs .745.
     - Next (on norm): seeds 3-5 per arm, mean +- SE closed / holo F1, plus a pose-share gate statistic at qtoyflex, qtoy2,
       IgM and r131.
     - Then a recommendation: keep on, gate on measured pose share (2 SE), or default off. Nothing reverted yet.
     - Gate statistic: the root latent's shift-error R^2 is .30-.35 on every IgM root (.002 with the profile), qtoy2 .034,
       qtoyflex .003, r131 0.
       - A significance gate fails: R^2 .003 clears 2 SE at 60k particles, and it needs true poses.
       - The truth-free share of latent evidence the profile removes separates by an order of magnitude: 75-87 % on IgM, <= 0
         elsewhere. That is the candidate gate.
   The native port is in progress (user: port first); then the r131 / 692k global-coordinate test and the end-to-end chain.
9. The r131 closed-corner gate: its readout passes on own-data imprint (corner leaf 1.03 / .54 own vs .93 / .59 on a fission
   re-read; random 5 % tails .87-.90). Options: keep the gate as is (it still separates arms consistently), or re-read every
   gate leaf on a fission copy with a threshold calibrated on honest readouts (e.g. against random tails of matched N).
   - DECIDED 2026-11-08 (user: "Yes, honest re-read"): the gate is read on a held-out fission copy, with the pass bar set from
     random tails of the same size. extremes-repro builds it (7.10 thread 8). Until it lands, own-data r131 corner readouts are
     provisional.

## 9. Standing rules and traps

Science:
- HARD GATE: flex works only at frozen poses; no latent recalculation or Woodbury algebra inside refinement unless rigorously
  proven useful.
- HARD GATE: the RNase P closed state must be retained (the r131 corner, the 692k closed leaves).
- No predefined region or multibody basis for RNase P; the tree is input-free (no seeding from prior knowledge); splits are data
  derived and the landscape must reproduce; science before speed; numerically close is enough, never chase bitwise identity.
- A new coupling must beat section 6's numbers on qtoy2, r131, qtoynull and olenull first.

- TRAP (2026-11-07): the RNase P denoised parts (`/mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part*.mrcs`,
  every r131 and 692k fixture) are ctf=flip, and the Python framework `$L/latent_states_v4.py` multiplies every stack by the
  SIGNED CTF (`X *= ctf`), so their signal past the first zero (~18 A at 1.7 um) enters inverted; `scale_noise`
  clips the resulting negative per-ring scales (rings past ~17 A) to zero weight, so in effect the band is cut at the first zero
  instead of flex's 15 A (rings 14-16 of 16 lost; high-defocus particles keep some inverted rings below the cut). Affected: every real-data
  run of the seam-free coordinate and the site dictionary (caches `$L/cache/f2r_shared_{r131,full692}.npy`); suspect: the other
  11-06 r131 caches (`f2r_fast_r131_r1_*`, `fp_fast_r131_r1`, `sp_mm2_r131_*`, `r131_r1`, `f2r_mm2_r131_seed`). Toys, olenull,
  IgM and inode dumps are ctf=yes. Use |CTF| for flipped stacks (memory rnasep-denoised-stacks-are-phase-flipped); the f2r
  session reran the r131 seam-free readout and the site dictionary with |CTF| (`$L/absctf/`): both move by <= .004, because their
  basis is band-limited at 15 A, so those verdicts stand. The trap matters wherever data past ~15 A is used (the klem boards at
  20-10 A changed with |CTF|). SIMPLE's own flex_pca coordinates and the r131 per-particle state-evidence gap (refine3D's
  log-evidence file) are unaffected. Fix `latent_states_v4.py` itself (|CTF| when the stk segment says flip) once no job runs it.

Cluster and data:
- Never delete cluster data (read, or `mv` aside); never signal SIMPLE processes by name; `scancel` by job id only; never rebuild or
  overwrite scripts under running jobs.
- `simple_exec` keys must follow `prg=`; a failed `chdir` does not abort (the run lands in the project root); finished trees prune
  split inputs.
- Memory sizing from `sacct` MaxRSS; wall time sized to the job.
- Partition (user 2026-11-08: "just use norm"): new jobs go to `-p norm`, submitted there directly so their QOS is norm, with
  `-c` <= 36 and `qsys_partition=norm` in every SIMPLE call. csbdevel stays usable too (per-user caps ~642 csbdevel, 600 norm).
  Running jobs are never moved.

Git:
- Push only to the fork; commits carry no Co-Authored-By trailer; build and fast-gate on batch2 before pushing code.
- Side sessions start with `use_worktree: false`: app-made worktrees fail on this repo ("shadow-unavailable", 4 failed starts on
  2026-11-07, and a retry reuses the failing setting). A code-writing session first runs `git -C /Users/atifao/SIMPLE worktree add
  --detach /Users/atifao/SIMPLE_<task> <commit>` and works only there; nothing is edited in `/Users/atifao/SIMPLE` itself.
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
