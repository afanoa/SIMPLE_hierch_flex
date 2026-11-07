#!/bin/bash
#SBATCH -J land_xdbasins
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 05:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_xdbasins.%j.out
# deconvolved latent basins (XD mixture through the flex noise model, modes and maximin barriers); the arguments go to land_xdbasins.py unchanged
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "[$(date '+%F %T') xdbasins] $*"
$PY -u $T/land_xdbasins.py "$@" 2>&1 | grep --line-buffered -v Warning
echo "[$(date '+%F %T') xdbasins] DONE"
