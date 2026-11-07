#!/bin/bash
#SBATCH -J fastgate_se
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/fast_gate.out
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:/mnt/nasapps/production/python/3.12/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/src && scripts/run_fast_gate.sh /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/build > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/fast_gate.log 2>&1; echo $? > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/fast_gate.exit
