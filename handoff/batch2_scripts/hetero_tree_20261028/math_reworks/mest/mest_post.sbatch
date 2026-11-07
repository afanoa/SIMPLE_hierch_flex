#!/bin/bash
#SBATCH -J mest_post
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest/logs/%x_%j.out
# stage-B posterior sharpness / calibration readout of the mu gate trees (text files only)
set -uo pipefail
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2
mkdir -p $M/runs/post
$PY -I $M/mest_post.py $M/runs/post/post.txt; echo "exit $?"
