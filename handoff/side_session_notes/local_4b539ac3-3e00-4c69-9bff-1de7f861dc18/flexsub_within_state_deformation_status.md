# Within-state deformation removal before reconstruction (flexsub): status

Owner: the flexsub side session (worktree laughing-herschel-157990). Cluster dir:
`/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/flexsub/` (scripts flexsub.py, flexsub_replay.sbatch,
flexsub_py.sbatch, flexsub_cut.py, flexsub_fsc.py; runs/py_<fixture>/ = correction diagnostics + flexsub.npz; runs/<key>_<arm>/ =
replays; stk/<fixture>/<arm>/all.mrcs = corrected stacks, stk/<fixture>/orig = symlink to the original stack for the control).

Method: at the node's stage-B poses, the latent model of latent_states_v4 (flex mean + 10 flex modes cut to 15 A, site_em scale/noise,
per-particle amplitude) with 2x-padded trilinear slices + gridding correction (the prototypes' nearest-neighbour cells are too crude
to subtract with). Class Gaussians of the stage-B labels fitted on one half, evaluated on the other; dz_i = E[z_i|y_i,c_i] - xi_c.
Subtract alpha_i beta CTF P_i(W dz_i) (alpha refitted per particle with the full latent model: alpha/a median 1.004, corr 0.999),
padded to the box, shifted back to the raw frame. Exactness: re-read planes = X - D to 5e-8, untouched beyond the flex band.
Arms: full (dz) and perp (dz with the between-state discriminant directions projected out in the Omega^-1 metric).
Replay: node recB input project -> update_project ptcl_old_root/ptcl_new_root -> saved recB, stage B, split calls (merge_master2,
SIMPLE_SEED 1, SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes).

## Correction diagnostics (runs/py_<fixture>/flexsub.log; flex band rings 2..16, noise-whitened, per particle)

fixture | |D|^2 / band residual | in-sample drop | held-out drop (other parity) / in-sample on same cells | truth R^2(s1 | dz) within states
---|---|---|---|---
qtoyflex root (3 states) | 0.02 / 851 (0.002 %) | 0.40 nats | 0.004 / 0.205 (ratio 0.02; 51 % of ptcls drop) | 0.019 / 0.023 / 0.057 (s1 sd 2.8-3.6); posterior mean over all arm ptcls 0.33
qtoy2 root (5 states) | 0.24 / 838 (0.03 %) | 1.56 nats | 0.097 / 0.780 (ratio 0.12; 59 %) | 0.010 / 0.169 / 0.060 / 0.151 (the merged-conformer states 2, 5); all arm 0.82

## Replays (runs/<key>_<arm>/readout/)

qtoyflex root: orig (uncorrected control) = production bit-for-bit (same label 1.000, split K* 3 27380/21183/11419).
arm | same label vs prod | closed F1 | holo F1 | kind NMI (split) | split
---|---|---|---|---|---
production | - | 0.169 (rec .727 pur .096) | 0.619 | 0.225 | K* 3: 27380/21183/11419
full | 0.991 | 0.166 (rec .715 pur .094) | 0.619 | 0.224 | K* 3: 27350/21212/11420
perp | 0.992 | 0.168 (rec .721 pur .095) | 0.619 | 0.225 | K* 3: 27371/21175/11436
Arm cut identical in all arms: s1 [-5.1, 3.4] / [-8.6, -1.6] / holo group [-7.1, 2.9] (p10-p90).

qtoy2 root: orig = production bit-for-bit (same label 1.000; split K* 5 21951/11451/10537/9413/6648).
arm | same label vs prod | truth NMI | closed F1 | holo F1 | split
---|---|---|---|---|---
production | - | 0.656 | 0.55 | 0.94 | K* 5
full | 0.997 | 0.656 | 0.55 | 0.94 | K* 5: 21928/11436/10564/9418/6654 (choose-rule margin 901 -> 756)

qtoynull correction: 2 stage-B states (latent Mahalanobis 5.4 between the arbitrary halves), |D|^2 0.01 / 851, in-sample 0.34 nats,
held-out 0.001 / 0.178 (ratio 0.00).
