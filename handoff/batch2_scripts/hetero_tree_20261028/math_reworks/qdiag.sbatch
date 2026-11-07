#!/bin/bash
#SBATCH -J mrw_proto
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
/mnt/nasapps/production/python/3.12/bin/python3 -I qdiag.py > qdiag_out.txt 2>&1
