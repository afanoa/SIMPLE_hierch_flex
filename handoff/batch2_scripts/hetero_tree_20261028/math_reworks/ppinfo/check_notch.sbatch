#!/bin/bash
#SBATCH -J ppinfo_chk
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo/logs/%x_%j.out
set -uo pipefail
P=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo; PY=/mnt/nasapps/production/python/3.12/bin/python3; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
N=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim2/runs/new_s1/P_L065/1_hetero_tree/node_r
echo "=== SIMPLE mean half maps (box 160) at 48"; $PY $P/check_notch.py $N/mean/recvol_state01_even.mrc $N/mean/recvol_state01_odd.mrc 48 $P/cache/L065s1_B88.npy
echo "=== recB state01 half maps (box 88) at 48"; $PY $P/check_notch.py $N/recB/recvol_state01_even.mrc $N/recB/recvol_state01_odd.mrc 48
echo "=== shift sign flipped: geomcheck + one posterior"; mkdir -p $P/runs/chk_shiftsign && cd $P/runs/chk_shiftsign
$PY $P/ppinfo_site.py --node $N --dump $P/dumps/L065s1 --site /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/lever/sites/c0_L065s1_site5 --mask $N/sites/site_05.mrc --dirs 1,3 --mskdiam 225.456 --out . --cache $P/cache/L065s1_B88.npy --cachebox 88 --nproc 1 --truth /mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim/fix/P_L065/truth.txt --rare 18 --labels file:/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/lever/sites/c0_L065s1_site5/site_post_z1lo.txt --geomcheck $N/mean/recvol_state01.mrc --shiftsign -1 2>&1 | grep -E "geomcheck|^\[.*\]  +[0-9]+ |N\* |summary|z1lo"
