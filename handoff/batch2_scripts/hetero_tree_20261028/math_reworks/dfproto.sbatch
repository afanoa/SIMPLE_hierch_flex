#!/bin/bash
#SBATCH -J mrw_dfproto
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
for nd in gates/base_r131/1_hetero_tree/node_r1 gates/base_r131_s1/1_hetero_tree/node_r1 gates/sub4_r131_s1/1_hetero_tree/node_r1 gates/sub4_r131/1_hetero_tree/node_r2 gates/sub4_qtoy2_s1/1_hetero_tree/node_r3 gates/sub4_qtoy2_s1/1_hetero_tree/node_r5 gates/base_qtoynull/1_hetero_tree/node_r gates/base_olenull/1_hetero_tree/node_r gates/base_micnull/1_hetero_tree/node_r; do
  /mnt/nasapps/production/python/3.12/bin/python3 -I drawfree_proto.py $nd 2>&1
done > dfproto_out.txt
