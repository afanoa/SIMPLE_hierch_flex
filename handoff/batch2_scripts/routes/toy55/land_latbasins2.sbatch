#!/bin/bash
#SBATCH -J land_latbasins2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_latbasins2.%j.out
# latent-density basins (ToMATo on the kNN density of flex latents); the arguments go to land_latbasins2.py unchanged
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
echo "[$(date '+%F %T') latbasins] $*"
$PY -u $T/land_latbasins2.py "$@" 2>&1 | grep --line-buffered -v Warning
echo "[$(date '+%F %T') latbasins] DONE"
