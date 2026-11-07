#!/bin/bash
#SBATCH -J mrw_draws
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%A_%a.out
# standalone state_split replays (instrumented build: production test logged, draw-free single pass decides), one subsample seed each
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks; B=$M/build_df/install
NODES=(base_r131:r1:232:c1:460 base_r131_s1:r1:232:c1:460 sub4_r131_s1:r1:232:c1:460 sub4_qtoy2_s1:r3:232:c1:132 sub4_qtoy2_s1:r5:232:c1:132 base_qtoynull:r:232:c1:132)
i=$((SLURM_ARRAY_TASK_ID % 6)); seed=$((SLURM_ARRAY_TASK_ID / 6 + 1))
IFS=: read -r t nd md pg nt <<< "${NODES[$i]}"
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=16 SIMPLE_SPLIT_CONTRAST=maps SIMPLE_SEED=1 SIMPLE_SPLIT_SEED=$seed
O=$M/draws/${t}_${nd}_seed$seed; mkdir -p $O; cd $O
cp $M/gates/$t/1_hetero_tree/node_$nd/split/s.simple s.simple
$B/bin/simple_exec prg=state_split projfile=s.simple mskdiam=$md pcg_mskfile=$M/gates/$t/1_hetero_tree/node_$nd/mask.mrc pgrp=$pg split_ntests=$nt nthr=16 mkdir=no > split.log 2>&1
echo "exit $?" >> split.log
