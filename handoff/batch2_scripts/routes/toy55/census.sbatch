#!/bin/bash
#SBATCH -J klem_census
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# fixed-point census of a two-class site EM on a toy with truth (census.py); env ARGS (census.py arguments)
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
/mnt/nasapps/production/python/3.12/bin/python3 -u /mnt/beegfs/elmlund/afan/routes/toy55/census.py ${ARGS:?} 2>&1 | grep -v Warning
