#!/bin/bash
#SBATCH -J pa_tail
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/logs/%x_%j.out
set -uo pipefail
A=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor
G=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
while read -r KEY NODE ARM SEED TRUTH XA XB; do
  [ "$ARM" = def ] || continue
  $PY $A/pa_tail.py $KEY $NODE $TRUTH $A/drift/${KEY}_n.txt $A/drift/${KEY}_b.txt
done < <(awk '!s[$1]++' $A/pa_list1.txt | awk '{$3="def"; print}')
