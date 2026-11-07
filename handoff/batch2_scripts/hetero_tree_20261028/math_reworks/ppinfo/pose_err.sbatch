#!/bin/bash
#SBATCH -J ppinfo_pose
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo/logs/%x_%j.out
P=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo; PY=/mnt/nasapps/production/python/3.12/bin/python3; export OMP_NUM_THREADS=1
for k in L065s1 L065s2 L08s1; do echo "=== $k"; f=$([ $k = L08s1 ] && echo P_L08 || echo P_L065); $PY $P/pose_err.py $P/dumps/$k /mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim/pool2/all.txt /mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim/fix/$f/truth.txt 18; done
