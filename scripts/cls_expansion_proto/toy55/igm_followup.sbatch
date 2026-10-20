#!/bin/bash
#SBATCH -J igm_followup
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Follow-ups on an IgM-like set once run_igm_sim has refined its poses (env MODE armsigm1): per-arm sector masks ->
# per-image likelihood ceiling of every arm at the TRUE poses (lp 20 and 10) and at the REFINED poses (lp 20) -> the
# pooled 2D route (abinitio2D + cls_expansion on the 5 best classes, synth_2droute.sbatch inline) read per arm.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:-armsigm1}; MK=$W/maps160/synth_$MODE/masks
log(){ echo "[$(date '+%F %T') igm_followup $MODE] $*"; }
[ -s $W/lat_$MODE/ptcl3D_r0.txt ] || { log "no refined poses yet ($W/lat_$MODE/ptcl3D_r0.txt)"; exit 1; }
[ -s $MK/all_arms.mrc ] || $PY $W/igm_arm_masks.py $MK || exit 1
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
for LP in 20 10; do echo "== ceiling at the TRUE poses, lp $LP"; $PY $W/arm_llr.py $W/sim$MODE $W/truth_$MODE.txt $MK/all_arms.mrc $W/maps160/synth_$MODE $W/armllr_$MODE/true_lp$LP --nsub 3000 --lp $LP --far 30 2>&1 | grep -v Warning | tail -14; done
echo "== ceiling at the REFINED poses, lp 20"; $PY $W/arm_llr.py $W/sim$MODE $W/truth_$MODE.txt $MK/all_arms.mrc $W/maps160/synth_$MODE $W/armllr_$MODE/refined_lp20 --nsub 3000 --lp 20 --far 30 --oris $W/lat_$MODE/ptcl3D_r0.txt 2>&1 | grep -v Warning | tail -14
log "pooled 2D route"; MODE=$MODE MSK=${MSK:-220} bash $W/synth_2droute.sbatch 2>&1 | tail -30
log DONE
