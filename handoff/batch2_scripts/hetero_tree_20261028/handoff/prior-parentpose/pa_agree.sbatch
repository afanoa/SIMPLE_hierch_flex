#!/bin/bash
#SBATCH -J pa_agree
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/logs/%x_%j.out
A=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
awk '!s[$1]++ {print $1, $2, $5}' ${LIST:-$A/pa_list1.txt} | while read -r KEY NODE TRUTH; do
  ls -d $A/replays/${KEY}_* > /dev/null 2>&1 || continue
  $PY $A/pa_agree.py $KEY $NODE $TRUTH
done
