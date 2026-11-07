#!/bin/bash
#SBATCH -J mrw_cmp
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
/mnt/nasapps/production/python/3.12/bin/python3 -I compare_sites.py dump_sites sub_sites > compare_sites.txt 2>&1
