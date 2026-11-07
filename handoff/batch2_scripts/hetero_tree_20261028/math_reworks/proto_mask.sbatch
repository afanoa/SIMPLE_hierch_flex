#!/bin/bash
#SBATCH -J mrw_proto
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
export OMP_NUM_THREADS=4
for nd in gates/base_r131/1_hetero_tree/node_r; do
  /mnt/nasapps/production/python/3.12/bin/python3 -I mask_disc_proto.py $nd --maps cheap 2>&1
done > proto_mask_out.txt
