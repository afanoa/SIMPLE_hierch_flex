#!/bin/bash
#SBATCH -J mrw_latent
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks; PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1
OUT=$M/gates/cb48_latent.txt; : > $OUT
for n in node_r node_r1 node_r2; do $PY -I $M/compare_latent.py $M/gates/bc_qtoy2 $M/gates/cb48_qtoy2 $n >> $OUT 2>&1; done
for n in node_r node_r1; do $PY -I $M/compare_latent.py $M/gates/cb_r131 $M/gates/cb48_r131 $n >> $OUT 2>&1; done
for n in node_r; do $PY -I $M/compare_latent.py $M/gates/head_qtoy2 $M/gates/bc_qtoy2 $n >> $OUT 2>&1; done
echo done >> $OUT
