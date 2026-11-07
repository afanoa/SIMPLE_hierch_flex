# prior-flexcorr: within-state deformation removal before reconstruction — VERDICT

**NO GAIN. Not worth porting.** Subtracting each particle's within-state continuous deformation (the node's own flex model, frozen
stage-B poses, computed once outside refinement) before recB -> stage B -> split changes nothing beyond seed noise on the three gate
roots: labels move 0.4-1.1 % (a seed change alone moves 2.9-3.9 %), truth F1 / NMI, the split and the state-map FSCs are unchanged,
and the null stays one leaf. Status: finished 2026-11-07 00:20 EST; no jobs open; no Fortran written.

## What was tested

- Correction (`flexsub.py`): at the node's stage-B poses, latent_states_v4's model (flex mean `node_r/mean/recvol_state01.mrc` + the
  node's 10 modes `fflex/flex_pca_pc0NN.mrc` cut to 15 A, site_em per-shell scale/noise, per-particle amplitude), with 2x-padded
  trilinear slices + gridding correction (the prototypes' nearest-neighbour cells are too crude to subtract with). Class Gaussians
  (tied Omega) of the production stage-B labels fitted on one half, evaluated on the other: dz_i = E[z_i | y_i, c_i] - xi_{c_i}.
  Subtracted alpha_i beta CTF_i P_i(W dz_i), alpha_i refitted per particle with the full latent model (alpha/a median 1.004, corr
  0.999), zero-padded to box 160, shifted back to the raw frame. Exactness: re-read planes = X - D to 5e-8; untouched beyond 15 A.
- Arms: `full` (dz) and `perp` (dz with the between-state discriminant directions projected out in the Omega^-1 metric, so the
  particle's latent state evidence is untouched). `orig` = uncorrected control through the same path (stack root = symlink to the
  original stack); `origs2` = the control at SIMPLE_SEED 2 (seed-noise scale).
- Replay (`flexsub_replay.sbatch`): copy of `node_r/recB/b.simple` (stage A's poses + folded labels) -> `simple_exec
  prg=update_project ptcl_old_root=<orig stack dir> ptcl_new_root=<arm stack dir>` -> the node's saved recB, stage B and split calls
  (merge_master2 build, SIMPLE_SEED 1, SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes) -> readouts. The orig control reproduces
  production bit for bit on all three fixtures (labels, split sizes, split statistic, map FSCs).

## Numbers

Correction diagnostics (flex band rings 2..16 = 117-14.6 A, noise-whitened sums per particle):

| fixture | subtracted power / band residual | in-sample drop | held-out drop (dz from the other checkerboard parity) vs in-sample on the same cells | truth R^2 of arm s1 on dz, within states |
|---|---|---|---|---|
| qtoyflex root (K 3) | 0.02 / 851 (0.002 %) | 0.40 nats | 0.004 vs 0.205 (ratio 0.02) | 0.019 / 0.023 / 0.057 (posterior mean over all arm ptcls: 0.33) |
| qtoy2 root (K 5) | 0.24 / 838 (0.03 %) | 1.56 nats | 0.097 vs 0.780 (ratio 0.12) | 0.010 / 0.169 / 0.060 / 0.151 (states 2, 5 merge conformers) |
| qtoynull root (K 2) | 0.01 / 851 (0.001 %) | 0.34 nats | 0.001 vs 0.178 (ratio 0.00) | - |

qtoyflex root (truth `routes/toy55/simqtoyflex/truth_qtoy.txt`):

| arm | same label vs production | closed F1 (recall, purity) | holo F1 | kind NMI | split | split z (crit 3.38) | mean FSC 15-10 A, states 1/2/3 |
|---|---|---|---|---|---|---|---|
| production = orig | 1.000 | 0.169 (.727, .096) | 0.619 | 0.225 | K* 3: 27380/21183/11419 | 27.1 | .8675/.8908/.8234 |
| origs2 | 0.971 | 0.168 | 0.620 | 0.226 | K* 3: 27372/21202/11408 | 28.5 | .8658/.8908/.8209 |
| full | 0.991 | 0.166 (.715, .094) | 0.619 | 0.224 | K* 3: 27350/21212/11420 | 26.5 | .8670/.8914/.8224 |
| perp | 0.992 | 0.168 (.721, .095) | 0.619 | 0.225 | K* 3: 27371/21175/11436 | 24.8 | .8673/.8909/.8234 |

How the arm is cut (identical in every arm, s1 p10-p90): child A [-5.1, 3.4] (holds closed, recall .72), child B [-8.6, -1.6], and a
holo child (47 % holo, 35 % arm particles over [-7.1, 2.9]): two overlapping cuts of the continuum.

qtoy2 root: NMI vs truth production/orig 0.656, full 0.656, perp 0.655; closed F1 0.55 and holo F1 0.94 in all; K* 5 in all
(full 21928/11436/10564/9418/6654 vs 21951/11451/10537/9413/6648); same label 0.997 / 0.996; split z 41.8 / 42.4 / 46.3; state-map
FSCs equal to +-0.0003.

qtoynull root: one child in every arm. Split z (K* before the guard): production = orig 4.68 (2, guard -> 1), origs2 3.43 (2 -> 1),
full 2.15 (1), perp 2.14 (1); same label 0.989 (full, perp), 0.961 (origs2). No anchoring of the null.

## Why

1. Per particle the within-state deformation is invisible: by the node's own model it is ~0.05 % of a particle's flex-band power, so
   the posterior mean is ~96 % shrunk to the state mean; the held-out gain (0.004 nats) matches that level. The flex latent sees the
   arm continuum only through the state labels (R^2 0.33 across states, 0.02-0.06 within).
2. In a linear latent model the within-state term is zero-mean around the state mean: subtracting it cannot un-blur a state map, only
   remove variance, bounded by that 0.05 %. Sharpening a continuum needs a nonlinear warp, which the frozen-pose gate excludes.
3. Where dz carries signal (qtoy2 states that merge conformers: held-out ratio 0.12, R^2 0.15-0.17) it is discrete sub-state
   evidence; subtracting it before a child's refinement would erase what the child needs to split. Neutral at the node, harmful in
   the recursion.

## What ran where (batch2)

Root: `/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/flexsub/` (scripts, lists, `logs/`, `runs/`, `stk/`, `cache/`).

- Corrections (`flexsub_py.sbatch`, LIST `list_py1.txt`): smoke (3000 ptcls) 55010347_1, qtoyflex 55010516_2, qtoy2 55010589_3,
  qtoynull 55011573_4 -> `runs/py_<tag>/` (flexsub.log, flexsub.npz with dz, dz_perp, dzA, posterior means, alpha, residual
  terms; summary.json).
- Corrected stacks: `stk/<fixture>/{full,perp}/all.mrcs` for qtoyflex, qtoy2, qtoynull (+ `stk/smoke_qtoyflex/`), ~6 GB each,
  44 GB in all; `stk/<fixture>/orig/all.mrcs` = symlinks to `/mnt/beegfs/elmlund/afan/routes/toy55/sim<fixture>/all.mrcs`.
  Nothing else reads them: move aside if space matters (never rm).
- Replays (`flexsub_replay.sbatch`): qtoyflex orig 55010348, full 55011413, perp 55011572 (`list_rep1.txt`); qtoy2 orig 55011887,
  full 55012019, perp 55012543, qtoynull full 55012624, orig + perp 55013182 (`list_rep2.txt`); seed-2 controls 55013511
  (`list_rep3.txt`) -> `runs/mm2_<fixture>_<arm>/` (recB, stageB with b_ptcl3D.txt and maps, split with split_stats.txt, readout/).
- Final readout (`flexsub_read.sbatch`): 55013760 -> `logs/fsub_read_55013760.out` (state-map FSCs of every arm, the arm-cut table).
  Note: `runs/mm2_qtoyflex_orig/readout/cut.txt` was written before a parser fix and its stage-B tables are wrong; the read log is
  authoritative.

## How to resume / reproduce

All from the batch2 login node, sbatch only (each replay is an 8-core master that submits 4 x 8-core parts; corrections 16 cores).

```bash
F=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/flexsub
# correction for one fixture (LIST line: TAG KEY NODE TRUTH MAXN); writes runs/py_<TAG>/ and stk/<TAG>/{full,perp}/
cd $F && sbatch --array=2 --export=ALL,LIST=$F/list_py1.txt flexsub_py.sbatch
# replay one arm (LIST line: KEY NODE ARM SEED OLDROOT NEWROOT TRUTH KINDTRUTH); refuses an existing runs/<KEY>_<ARM>
cd $F && sbatch --array=2 --export=ALL,LIST=$F/list_rep1.txt flexsub_replay.sbatch
# FSC + arm-cut readout over every arm present
cd $F && sbatch flexsub_read.sbatch
```

A new fixture needs its dumps under `math_reworks/latent/dumps/<KEY>/` (n_ptcl3D.txt, b_ptcl3D.txt, n_stk.txt), a TAG line in a py
list, an `stk/<TAG>/orig` symlink for the control, and replay lines pointing NEWROOT at `stk/<TAG>/<arm>`. The uncorrected replay from
recB reproduces production bit for bit, a cheap harness for any other per-particle image edit.

Files here: flexsub.py, flexsub_cut.py, flexsub_fsc.py, flexsub_py.sbatch, flexsub_replay.sbatch, flexsub_read.sbatch, the four
lists (identical to the copies that ran), NOTES.md (the side-session notes).
