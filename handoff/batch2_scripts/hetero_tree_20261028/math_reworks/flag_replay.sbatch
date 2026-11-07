#!/bin/bash
#SBATCH -J mrw_flag
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%A_%a.out
# standalone state_split replays (the flag before the guard): TASK = tree:node:mskdiam:pgrp:ntests
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks; B=$M/build_r2i/install
TASKS=(base_qtoynull:r:232:c1:132 base_r131:r:232:c1:460 base_qtoy2:r:232:c1:132 base_qtoy2:r1:232:c1:132 base_olenull:r:240:c2:52 base_micnull:r:200:c1:36)
IFS=: read -r t nd md pg nt <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=16 SIMPLE_SPLIT_CONTRAST=maps SIMPLE_SEED=1
O=$M/flag2/${t}_$nd; mkdir -p $O; cd $O
cp $M/gates/$t/1_hetero_tree/node_$nd/split/s.simple s.simple
$B/bin/simple_exec prg=state_split projfile=s.simple mskdiam=$md pcg_mskfile=$M/gates/$t/1_hetero_tree/node_$nd/mask.mrc pgrp=$pg split_ntests=$nt nthr=16 mkdir=no > split.log 2>&1
echo "exit $?" >> split.log
