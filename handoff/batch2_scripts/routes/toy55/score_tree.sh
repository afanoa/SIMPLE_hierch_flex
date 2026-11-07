#!/bin/bash
# score_tree.sh RUN_DIR KIND [PAIR_RUN_DIR]: submits score_tree.sbatch (norm, 32 cpus) for the one-page usefulness scorecard of a
# hetero_tree run; KIND rnasep | cdisc | igm | cahra | other. The page lands in rnasep_klem/scorecard/<run>/scorecard.txt.
R=${1:?usage: score_tree.sh RUN_DIR KIND [PAIR_RUN_DIR]}; KD=${2:?KIND}; PR=${3:-}
TREE=$(cd $R && pwd) KIND=$KD PAIR=$PR sbatch -J score_$(basename $R) /mnt/beegfs/elmlund/afan/routes/toy55/score_tree.sbatch
