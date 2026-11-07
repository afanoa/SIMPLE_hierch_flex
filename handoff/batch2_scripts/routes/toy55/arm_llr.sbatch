#!/bin/bash
#SBATCH -J synth_armllr
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# per-image likelihood ceiling of the hinged-arm set at known poses (env MODE, REG = region mask, NSUB 3000, LP 10)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:?}; REG=${REG:?}; O=$W/armllr_$MODE/$(basename $REG .mrc); [ -n "${ORIS:-}" ] && O=${O}_refined; [ "${LP:-10}" != 10 ] && O=${O}_lp${LP}; mkdir -p $O
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
echo "== $MODE / $(basename $REG): per-image likelihood ceiling at the true poses"
TR=$W/truth_$MODE.txt; [ -f $TR ] || TR=none; MAPS=$W/maps160/synth_$MODE; [ -d $MAPS ] || MAPS=${MAPS_DIR:-$W/maps160/docked_all}
$PY $W/arm_llr.py $W/sim$MODE $TR $REG $MAPS $O --nsub ${NSUB:-3000} --lp ${LP:-10} --far ${FAR:-30} ${ORIS:+--oris $ORIS} 2>&1 | grep -v Warning
echo "[$(date '+%F %T') armllr] DONE"
