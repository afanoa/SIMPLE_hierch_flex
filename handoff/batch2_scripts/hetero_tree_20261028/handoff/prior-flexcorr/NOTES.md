# Within-state deformation removal before reconstruction (flexsub): VERDICT NO GAIN, do not port

Owner: the flexsub side session (worktree laughing-herschel-157990). Cluster dir:
`/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/flexsub/` (flexsub.py = correction, flexsub_py.sbatch,
flexsub_replay.sbatch = replay from recB, flexsub_cut.py = qtoyflex arm-cut readout, flexsub_fsc.py + flexsub_read.sbatch = state-map
FSCs; lists list_py1 / list_rep1-3). runs/py_<fixture>/flexsub.log + flexsub.npz (dz, dz_perp, posterior means, alpha, residual drops);
runs/<key>_<arm>/ = replays with readout/; stk/<fixture>/<arm>/all.mrcs = corrected stacks (~6 GB each, 8 + 2 smoke; nothing else
reads them: move aside if space matters), stk/<fixture>/orig = symlink to the original stack (control).

Method: at the node's stage-B poses, latent_states_v4's model (flex mean + the node's 10 flex modes cut to 15 A, site_em scale/noise,
per-particle amplitude) with 2x-padded trilinear slices + gridding correction. Class Gaussians of the stage-B labels fitted on one half,
evaluated on the other; dz_i = E[z_i|y_i,c_i] - xi_c. Subtract alpha_i beta CTF P_i(W dz_i) (alpha refitted with the full latent model:
alpha/a median 1.004, corr 0.999), padded to the box, shifted to the raw frame; re-read planes = X - D to 5e-8, untouched beyond 15 A.
Arms: full (dz) and perp (dz with the between-state discriminant directions projected out, Omega^-1 metric).
Replay: the node's recB input project -> update_project ptcl_old_root/ptcl_new_root -> saved recB, stage B, split calls (merge_master2,
SIMPLE_SEED 1, SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes). The uncorrected control (orig) reproduces production bit for bit
on all three fixtures (labels, split sizes, split statistic, map FSCs); origs2 = the control at SIMPLE_SEED 2 (seed-noise scale).

## Correction diagnostics (flex band rings 2..16, noise-whitened sums per particle)

fixture | |D|^2 / band residual | in-sample drop | held-out drop (dz from the other checkerboard parity) vs in-sample on the same cells | truth R^2(s1 | dz) within states
---|---|---|---|---
qtoyflex root (K 3) | 0.02 / 851 (0.002 %) | 0.40 nats | 0.004 vs 0.205 (ratio 0.02; 51 % of ptcls drop) | 0.019 / 0.023 / 0.057 (s1 sd 2.8-3.6); posterior mean over all arm ptcls 0.33
qtoy2 root (K 5) | 0.24 / 838 (0.03 %) | 1.56 nats | 0.097 vs 0.780 (ratio 0.12; 59 %) | 0.010 / 0.169 / 0.060 / 0.151 (states 2, 5 merge conformers); all arm 0.82
qtoynull root (K 2) | 0.01 / 851 (0.001 %) | 0.34 nats | 0.001 vs 0.178 (ratio 0.00; 51 %) | -
Model-implied within-state deformation power per particle (from the in-sample drop, small-lambda limit): qtoyflex ~0.4 whitened units of
851 (0.05 %), so the posterior mean is ~96 % shrunk to the state mean; the held-out gain matches that tiny level.

## Replays

qtoyflex root | same label vs prod | closed F1 (rec, pur) | holo F1 | kind NMI | split K*, sizes | split z (crit 3.38) | state-map mean FSC 15-10 A (s1/s2/s3)
---|---|---|---|---|---|---|---
production = orig | 1.000 | 0.169 (.727, .096) | 0.619 | 0.225 | 3: 27380/21183/11419 | 27.1 | .8675/.8908/.8234
origs2 (seed 2) | 0.971 | 0.168 | 0.620 | 0.226 | 3: 27372/21202/11408 | 28.5 | .8658/.8908/.8209
full | 0.991 | 0.166 (.715, .094) | 0.619 | 0.224 | 3: 27350/21212/11420 | 26.5 | .8670/.8914/.8224
perp | 0.992 | 0.168 (.721, .095) | 0.619 | 0.225 | 3: 27371/21175/11436 | 24.8 | .8673/.8909/.8234
Arm cut identical in every arm (s1 p10-p90): child A [-5.1, 3.4] (holds closed, recall .72), child B [-8.6, -1.6], holo child (47 % holo,
35 % arm over [-7.1, 2.9]): two overlapping cuts of the continuum.

qtoy2 root | same label | truth NMI | closed F1 | holo F1 | split | z
---|---|---|---|---|---|---
production = orig | 1.000 | 0.656 | 0.55 | 0.94 | K* 5: 21951/11451/10537/9413/6648 | 41.8
full | 0.997 | 0.656 | 0.55 | 0.94 | K* 5: 21928/11436/10564/9418/6654 | 42.4
perp | 0.996 | 0.655 | 0.55 | 0.94 | K* 5: 21950/11435/10551/9409/6655 | 46.3
State-map FSCs equal to +-0.0003 in every band.

qtoynull root | same label | split | z (K* before the guard)
---|---|---|---
production = orig | 1.000 | 1 child | 4.68 (2, guard -> 1)
origs2 | 0.961 | 1 child | 3.43 (2, guard -> 1)
full | 0.989 | 1 child | 2.15 (1)
perp | 0.989 | 1 child | 2.14 (1)
No anchoring: the corrected nulls sit further below z_crit than production; the shift is within the null's seed spread.

## Reading
1. Per particle the within-state deformation is invisible: ~0.05 % of a particle's flex-band power by the node's own model, held-out
   gain 0.004 nats, within-state R^2 of the arm coordinate 0.02-0.06. The latent sees the arm only through the state labels (R^2 0.33).
2. In a linear latent model the within-state term is zero-mean around the state mean: subtracting it cannot un-blur a state map, only
   remove variance, bounded by that 0.05 %. Measured: map FSCs move <= 0.0006, below seed noise; labels move 0.4-1.1 % (seed 2: 2.9-3.9 %).
3. Where dz does carry signal (qtoy2 states that merge conformers: held-out ratio 0.12, R^2 0.15-0.17) it is discrete sub-state signal;
   subtracting it before a child's refinement would erase what the child needs to split. A port would be neutral at the node and
   harmful in the recursion. Sharpening a continuum needs a warp (nonlinear), not a subtraction, and that is outside the frozen-pose gate.
No Fortran written.
