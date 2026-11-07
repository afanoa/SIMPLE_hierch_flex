#!/bin/bash
#SBATCH -J shrinkBc
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB/logs/%x_%j.out
# Two shrinkB_run.sbatch passes in one allocation: MODE1 over LIST1, then MODE2 over LIST2 (one slot of the 100-core budget).
set -uo pipefail
: ${MODE1:?} ${LIST1:?} ${MODE2:?} ${LIST2:?}
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB
MODE=$MODE1 LIST=$LIST1 bash $M/shrinkB_run.sbatch
MODE=$MODE2 LIST=$LIST2 bash $M/shrinkB_run.sbatch
