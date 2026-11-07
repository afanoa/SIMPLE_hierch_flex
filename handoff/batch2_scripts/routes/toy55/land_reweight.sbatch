#!/bin/bash
#SBATCH -J land_reweight
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_reweight.%j.out
# ensemble reweighting of a scored path. env SC OUT ARGS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
$PY -u $T/land_reweight.py --scores ${SC:?} --out ${OUT:?} ${ARGS:-} 2>&1 | grep --line-buffered -v "Warning\|warn("
echo "[$(date '+%F %T')] land_reweight DONE"
