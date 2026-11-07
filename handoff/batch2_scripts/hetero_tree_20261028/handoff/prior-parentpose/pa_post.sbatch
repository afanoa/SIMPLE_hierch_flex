#!/bin/bash
#SBATCH -J pa_post
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/logs/%x_%j.out
# per replay: the full in -> B drift and the tail by truth class (the replay script's drift.txt keeps only the geodesic lines)
set -uo pipefail
A=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
while read -r KEY NODE ARM SEED TRUTH XA XB; do
  R=$A/replays/${KEY}_${ARM}_s$SEED
  [ -s $R/stageB/b_ptcl3D.txt ] || continue
  echo "######## ${KEY}_${ARM}_s$SEED"
  $PY $A/pa_drift.py ${KEY}_${ARM}_s$SEED $A/drift/${KEY%_r?}_r_b.txt $A/drift/${KEY}_n.txt $R/stageA/a_ptcl3D.txt $R/stageB/b_ptcl3D.txt 1.464 $NODE/rows.txt 2>&1 | sed -n '/^in->B/,/shift/p'
  $PY $A/pa_tail.py ${KEY}_${ARM}_s$SEED $NODE $TRUTH $A/drift/${KEY}_n.txt $R/stageB/b_ptcl3D.txt 2>&1 | sed -n '2,5p'
done < ${LIST:-$A/pa_list1.txt}
