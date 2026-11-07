#!/bin/bash
#SBATCH -J synth_llrtol
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# pose-error tolerance of the per-image likelihood ceiling (env MODE, REG, MAPS_DIR, LP 20, SIGS "2 5 10"): the refined
# table's error against the truth, then arm_llr.py at the truth perturbed by each sigma
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:?}; REG=${REG:?}; LP=${LP:-20}
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
TR=$W/truth_$MODE.txt; [ -f $TR ] || TR=none; MAPS=${MAPS_DIR:-$W/maps160/docked_all}
echo "== $MODE pose-error tolerance at lp $LP"
$PY $W/pose_perturb.py $W/sim$MODE/all.txt $W/lat_$MODE/ptcl3D_r0.txt $W/sim$MODE/gt_particles.txt $W/sim$MODE/perturbed ${SIGS:-2 5 10} 2>&1 | grep -v Warning
for s in ${SIGS:-2 5 10}; do O=$W/armllr_$MODE/$(basename $REG .mrc)_sig${s}_lp$LP; mkdir -p $O; echo "== sigma $s deg"
  $PY $W/arm_llr.py $W/sim$MODE $TR $REG $MAPS $O --nsub ${NSUB:-2000} --lp $LP --oris $W/sim$MODE/perturbed_sig$s.txt 2>&1 | grep -a "accuracy 0\|d. median\|pooled   10 \|arm 1: LLR"; done
echo "[$(date '+%F %T') llrtol] DONE"
