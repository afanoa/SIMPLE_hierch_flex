#!/bin/bash
#SBATCH -J land_morph
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_morph.%j.out
# optimal-transport morph between two maps. env MA MB OUT [ARGS]
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
$PY -u $T/land_morph.py ${MA:?} ${MB:?} ${OUT:?} ${ARGS:-} 2>&1 | grep --line-buffered -v "Warning\|warn("
echo "[$(date '+%F %T')] land_morph DONE"
