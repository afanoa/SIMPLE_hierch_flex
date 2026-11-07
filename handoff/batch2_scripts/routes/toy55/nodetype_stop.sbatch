#!/bin/bash
#SBATCH -J nodetype_stop
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Stop test with sector-split noise (nodetype_stop.py) on every cheap-children dir of the dev set and, with SET=all, the real nodes.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
for d in null_50k_rc:root_null_50k_rc null_50k_r:root_null_50k_r two_50k_r:root_two_50k_r two_50k_rc:root_two_50k_rc s30a_r11:s30a_r11 s30a_r21:s30a_r21 s30a_r22:s30a_r22 s30a_r12:s30a_r12 s30a_r1:s30a_r1 s30a_r:s30a_r s30a_r2:s30a_r2; do
  $PY -u $T/nodetype_stop.py $K_/cheapkids_${d#*:} --tag "[${d%%:*}]" 2>&1 | grep --line-buffered -v Warning; done
if [ "${SET:-dev}" = all ]; then for u in r1 r21 r22 r2 r; do $PY -u $T/nodetype_stop.py $K_/cheapkids_f692rf_$u --mskdiam 232 --tag "[f692rf_$u]" 2>&1 | grep --line-buffered -v Warning; done; fi
echo "[$(date '+%F %T') nodetype_stop] DONE"
