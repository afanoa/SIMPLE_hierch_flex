#!/bin/bash
#SBATCH -J klem_leaf_bfac
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# sharpness vs structure of the leaf states (env RUN)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; RUN=${RUN:?}; echo "== leaf B-factors of $RUN"; $PY $W/leaf_bfac.py $RUN $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning; echo "[$(date +%T) bfac] DONE"
