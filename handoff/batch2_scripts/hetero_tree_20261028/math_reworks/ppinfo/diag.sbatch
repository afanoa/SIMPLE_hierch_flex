#!/bin/bash
#SBATCH -J ppinfo_diag
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=20G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo/logs/%x_%j.out
set -uo pipefail
P=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo; PY=/mnt/nasapps/production/python/3.12/bin/python3; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
N=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim2/runs/new_s1/P_L065/1_hetero_tree/node_r
$PY $P/diag_shells.py --node $N --dump $P/dumps/L065s1 --site /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/lever/sites/c0_L065s1_site5 --mask $N/sites/site_05.mrc --cache $P/cache/L065s1_B88.npy --vol $N/mean/recvol_state01.mrc 2>&1 | grep -v "^\[.*\] *z"
