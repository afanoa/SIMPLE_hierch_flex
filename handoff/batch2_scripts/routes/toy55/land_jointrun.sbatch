#!/bin/bash
#SBATCH -J land_joint
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_joint.%j.out
# joint map + particle fit on a synthetic occval run (truth basis). env RUN, PEN
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RUN=${RUN:?}
$PY -u $T/land_joint.py --run $RUN --kind e1b --sim $T/simswing30 --mr full --cache $K_/lsamp/caches/swing30 --poses $K_/lsamp/s30a/poses_combined.txt \
   --basis $K_/landscape/e0_s30/truth/leaves.txt --out $RUN/joint --lam_pen ${PEN:-1e-3} --nw 12 2>&1 | grep --line-buffered -v "Warning\|warn("
echo "[$(date '+%F %T')] land_joint $RUN DONE"
