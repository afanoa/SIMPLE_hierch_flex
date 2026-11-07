#!/bin/bash
#SBATCH -J nodemask_chk
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Focus-mask sources vs truth at every rflex_s30a node (nodemask_check.py): flex-state cross-half variance (current), cheap children's
# cross-half variance, flex PC maps; scored against the node's true real-space variance.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
RD=$K_/rflex_s30a
for u in ${NODES:-r r1 r2 r11 r12 r21 r22}; do ROWS=$(awk -v u=$u '$1 == u {print $4}' $RD/queue.txt)
  $PY $T/nodemask_check.py --node $RD/node_$u --kids $K_/cheapkids_s30a_$u --rows $ROWS --gt $T/simswing30/gt_particles.txt --tmaps $T/simswing30/maps --tag "[$u]" 2>&1 | grep -v Warning; done
echo "[$(date '+%F %T') nodemask_check] DONE"
