#!/bin/bash
#SBATCH -J land_e0b
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_e0b.%j.out
# smooth cross-half deconvolution of the oracle scores (true and node poses)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; O=$A/rnasep_klem/landscape/e0_s30
for v in ${VARS:-true node}; do $PY -u $T/land_e0b.py --scores $O/sc_$v --gt $T/simswing30/gt_particles.txt --out $O/b_$v 2>&1 | grep -v "Warning\|warn("; done
echo "[$(date '+%F %T')] land_e0b DONE"
