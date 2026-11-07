#!/bin/bash
#SBATCH -J land_basins
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_basins.%j.out
# basin trees: swing30 calibration (morph-path scores, lam 0.01 and 0.001), then real r1 on both paths
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; L=/mnt/beegfs/elmlund/afan/rnasep_klem/landscape
for lam in 0.01 0.001; do echo "== swing30 morph path, lam $lam"; $PY -u $T/land_basins.py --scores $L/ps_s30_morph.npz --gt $T/simswing30/gt_particles.txt --lam $lam 2>&1 | grep --line-buffered -v Warning; done
for p in 18_49 18_28; do for lam in 0.01 0.001; do echo "== real r1 path $p, lam $lam"; $PY -u $T/land_basins.py --scores $L/ps_r1_$p.npz --lam $lam 2>&1 | grep --line-buffered -v Warning; done; done
echo "[$(date '+%F %T')] land_basins DONE"
