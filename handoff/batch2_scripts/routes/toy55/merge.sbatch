#!/bin/bash
#SBATCH -J klem_merge
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}
$PY $W/merge_states.py $R $W/maps160/toy_consensus.mrc $R/merge_tree --lp ${LP:-10} 2>&1 | grep -v Warning; echo "MERGE DONE $RTAG"
