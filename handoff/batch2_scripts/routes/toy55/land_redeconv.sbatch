#!/bin/bash
#SBATCH -J land_redeconv
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_redeconv.%j.out
# smoothness-rule comparison of finished land_occ runs. env RUNS (space separated), KIND
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
for r in ${RUNS:?}; do echo "== $r"; [ -s $r/occ/occ.json ] || { echo "missing $r/occ/occ.json"; continue; }
  if [ ${KIND:-e1b} = e1b ]; then $PY -u $T/land_redeconv.py $r --truthocc $K_/landscape/occ_truth_table_s30.txt 2>&1 | grep --line-buffered -v Warning
  else $PY -u $T/land_redeconv.py $r --kind real 2>&1 | grep --line-buffered -v Warning; fi; done
echo "[$(date '+%F %T')] land_redeconv DONE"
