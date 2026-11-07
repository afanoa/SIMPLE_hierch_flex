#!/bin/bash
#SBATCH -J land_xdcompare
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 05:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_xdcompare.%j.out
# deconvolved latent landscape compared across runs (common frame, Bhattacharyya vs half ceiling, marginals); the arguments go to land_xdcompare.py unchanged
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "[$(date '+%F %T') xdcompare] $*"
$PY -u $T/land_xdcompare.py "$@" 2>&1 | grep --line-buffered -v Warning
echo "[$(date '+%F %T') xdcompare] DONE"
