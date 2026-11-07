#!/bin/bash
#SBATCH -J land_realpath
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_realpath.%j.out
# validated reweighting (lmarg, tau split, lam 0.01) + few-state test of one real scored path. env SC OUTB
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
$PY -u $T/land_reweight.py --scores ${SC:?} --out ${OUTB:?}_lmarg_split_0.01 --real --variant lmarg --tau split --lam 0.01 --nm 932496 --nboot 20 2>&1 | grep --line-buffered -v "Warning\|warn("
$PY -u $T/land_twostate.py --scores $SC 2>&1 | grep --line-buffered -v Warning
echo "[$(date '+%F %T')] land_realpath DONE"
