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

Each writes to `$H/handoff/<name>/` (VERDICT.md, or STATUS.md with job ids and resume commands). State at 17:50 EST:
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
  - Part 4b (queued): GAPw per block on the knob-free pursuit's p3 blocks, whose extra blocks Part 4 never saw.
  - Part 5 (running; design approved): a fission saturation test replacing .28. Fission in memory (tau 1) on the framework's
    cached Fourier planes; one y1 score orders the particles into four groups (each side of the input split halved at its
    median); y2 maps give cross-half squared steps q_A, q_mid, q_B. DISCRETE iff z(q_mid) >= 3 and one side has q_side / q_mid
    + 2 SE < 1: equal steps mean uniform occupancy along the path whatever its curvature, and a homogeneous state saturates.
    Expected r: closed side ~0, curved arm continuum ~1, IgM 10 % a1 ~.5. Known limit: a continuum peaked at one end reads as a
    state on its dense side (a subsampled qtoyflex control measures it). Runs on qtoyflex (closed vs rest, arm median split,
    a realistic closed split) and IgM 10 % a1 fix / null.
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
    (`src/main/flex/fit/simple_flex_pca_embed.f90`), default on, `SIMPLE_FLEX_POSE_PROFILE=no` turns it off. Build
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
  (eigenvalues up to 12 %). Neutral to slightly negative on the gates; a third qtoy2 seed pair is running. Not committed.
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
   (per-group spectrum; r131's stacks are non-stationary) before the gates; the band split is the fallback.
7. Pose profile in flex's embedding: done, `5f80720e0` on `hetero-tree-wip` after the gates passed (7.8). The flex-precision
   fix (stacked patch, 7.8): correct in principle and closer to truth at the root latent, but neutral to slightly negative on the
   gates (qtoy2 closed merged on 1 of 2 seeds); commit or hold after the third qtoy2 seed pair.
8. Site layer: sparse PPCA's pursuit (knob-free stop: 2 SE + solvent controls, leave-one-out, support components) + one masked
   round, the best data-only chain (7.8 igm-refs Part 3), with GAPw accepting each split (Part 4: 9/10 power, 0/10 false
   accepts; valid for IgM-like short swings, not for qtoyflex's curved arm continuum or its closed state next to it). The native
   port is in progress (user: port first); then the r131 / 692k global-coordinate test and the end-to-end chain.
9. The r131 closed-corner gate: its readout passes on own-data imprint (corner leaf 1.03 / .54 own vs .93 / .59 on a fission
   re-read; random 5 % tails .87-.90). Options: keep the gate as is (it still separates arms consistently), or re-read every
   gate leaf on a fission copy with a threshold calibrated on honest readouts (e.g. against random tails of matched N).

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
