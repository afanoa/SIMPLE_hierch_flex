#!/bin/bash
#SBATCH -J rnasep_umaps_docking
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# UMAPs for the user: (1) the family-1-wide focal board (409k) coloured by the docking rank sets (1 = within-10k closed
# set, 2 = rest of the docked top 30k, 3 = next 20k, 4 = undocked bottom 30k, 0 = the rest); (2) the within-docked board
# (30k) with the same colours; the usual density / z1 / z2 / junk / tilt / defocus panels alongside.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; F=$W/figs_now; mkdir -p $F; cd $W; export OMP_NUM_THREADS=16
$PY - <<PYEOF
import numpy as np
a = np.loadtxt('$KF/closed_focal/lever/states.txt', dtype=int); b = np.loadtxt('$KF/closed_focal_docked30k/rank/states.txt', dtype=int); lab = np.zeros(len(a), int)
lab[a == 3] = 4; lab[a == 2] = 3; lab[a == 1] = 2; lab[b == 1] = 1; np.savetxt('$F/docking_rank_labels.txt', lab, fmt='%d'); print('labels:', np.bincount(lab))
PYEOF
$PY $W/klem_umap.py $KF/closed_focal/board $P/meta_fb.npz --labels $F/docking_rank_labels.txt --label_title "docking rank sets (1 closed 10k, 2 rest of top 30k, 3 next 20k, 4 undocked)" --label_name "dock" --out $F/umap_real_focal_fam1_docking.png 2>&1 | grep -v Warning | tail -3
$PY $W/klem_umap.py $KF/closed_focal_docked30k/board $P/meta_fb.npz --labels $F/docking_rank_labels.txt --label_title "docking rank sets inside the docked 30k (1 = closed 10k)" --label_name "dock" --out $F/umap_real_within_docked30k.png 2>&1 | grep -v Warning | tail -3
echo "[$(date '+%F %T')] umaps DONE"
