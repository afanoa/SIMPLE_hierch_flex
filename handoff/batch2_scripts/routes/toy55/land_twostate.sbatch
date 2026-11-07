#!/bin/bash
#SBATCH -J land_twostate
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_twostate.%j.out
# two-state vs path test (synthetic morph-path scores, then real r1) and the reweighting tables
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; L=/mnt/beegfs/elmlund/afan/rnasep_klem/landscape

$PY -u $T/land_twostate.py --scores $L/ps_s30_morph.npz --gt $T/simswing30/gt_particles.txt 2>&1 | grep --line-buffered -v Warning
$PY -u $T/land_twostate.py --scores $L/ps_r1_18_49.npz 2>&1 | grep --line-buffered -v Warning
echo "[$(date '+%F %T')] land_twostate DONE"
