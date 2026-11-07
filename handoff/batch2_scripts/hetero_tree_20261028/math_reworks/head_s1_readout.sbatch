#!/bin/bash
#SBATCH -J mrw_s1read
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks; export PYTHONDONTWRITEBYTECODE=1
{ for a in head; do echo "## ${a}_r131_s1"; grep -E 'HETERO_TREE node [a-z0-9]*: (split group|states)' gates/${a}_r131_s1/tree.log | cut -c1-140; grep -h 'strongest kept pole' gates/${a}_r131_s1/tree.log | cut -c1-200; done
  /mnt/nasapps/production/python/3.12/bin/python3 -I /mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts/closed_nodes2.py gates/s1_readout df_s1=gates/head_r131_s1/1_hetero_tree 2>&1 | grep -v -E 'Deprecation|V = np'
} > gates/head_s1_readout.txt 2>&1
