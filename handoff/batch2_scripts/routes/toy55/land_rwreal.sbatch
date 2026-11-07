#!/bin/bash
#SBATCH -J land_rwreal
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_rwreal.%j.out
# real-data reweighting over a set of settings. env SC OUTB (output base), SETS ("variant:tau:lam ...")
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
for s in ${SETS:?}; do IFS=: read v tau lam <<< "$s"; echo "== $v tau $tau lam $lam"
  $PY -u $T/land_reweight.py --scores ${SC:?} --out ${OUTB:?}_${v}_${tau}_${lam} --real --variant $v --tau $tau --lam $lam --nm 932496 --nboot 20 2>&1 | grep --line-buffered -v "Warning\|warn("; done
echo "[$(date '+%F %T')] land_rwreal DONE"
