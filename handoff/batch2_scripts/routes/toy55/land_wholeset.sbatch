#!/bin/bash
#SBATCH -J land_wholeset
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_wholeset.%j.out
# whole-set (junk-free) occupancy along conformer 18 -> 49: validated reweighting, basin tree, few-state test
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; L=/mnt/beegfs/elmlund/afan/rnasep_klem/landscape; SC=$L/ps_all_18_49.npz; J=$L/rows_junk_f692st.txt
for lam in 0.01 0.001; do $PY -u $T/land_reweight.py --scores $SC --out $L/rwreal_all_18_49_lmarg_split_${lam}_nojunk --real --variant lmarg --tau split --lam $lam --nm 932496 --nboot 20 --exclude $J 2>&1 | grep --line-buffered -v "Warning\|warn("; done
for lam in 0.01 0.001; do echo "== basins lam $lam"; $PY -u $T/land_basins.py --scores $SC --lam $lam --thr 0.5,0.6,0.75 --exclude $J 2>&1 | grep --line-buffered -v Warning; done
$PY -u $T/land_twostate.py --scores $SC --exclude $J 2>&1 | grep --line-buffered -v Warning
echo "[$(date '+%F %T')] land_wholeset DONE"
