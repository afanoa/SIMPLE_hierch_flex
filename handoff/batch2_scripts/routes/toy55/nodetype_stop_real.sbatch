#!/bin/bash
#SBATCH -J nodetype_stopr
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Stop test with sector-split noise (nodetype_stop.py) on the real rflex_f692rf nodes (mskdiam 232).
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
for u in r r1 r2 r21 r22; do $PY -u $T/nodetype_stop.py $K_/cheapkids_f692rf_$u --mskdiam 232 --tag "[f692rf_$u]" 2>&1 | grep --line-buffered -v Warning; done
echo "[$(date '+%F %T') nodetype_stop real] DONE"
