#!/bin/bash
#SBATCH -J igm2_route2d
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The pooled 2D route on an IgM-like set in ONE job (truth only in readouts): synth_2droute.sbatch inline (abinitio2D, 100
# classes; cls_expansion sub-classes on the 5 best classes) -> sub-class labels printed with the stage1 build (rd2d.sbatch;
# the main build cannot read the expansion project's records) -> subclass_angle_readout.py per Fab (in-plane theta,
# out-of-plane phi), next to the same readout of the 5-arm set's sub-classes. env MODE (armsigm2), MSK (220), REF (armsigm1)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; export MODE=${MODE:-armsigm2} MSK=${MSK:-220}; REF=${REF:-armsigm1}
log(){ echo "[$(date '+%F %T') igm2_route2d $MODE] $*"; }
log "2D route (abinitio2D + cls_expansion)"; bash $W/synth_2droute.sbatch > $W/igm2_${MODE}_route2d.log 2>&1; grep -a -E "2D done|cls_expansion|FAILED" $W/igm2_${MODE}_route2d.log | tail -4
log "sub-class labels (stage1 build)"; bash $W/rd2d.sbatch > $W/igm2_${MODE}_rd2d.log 2>&1; grep -a -E "^dump|^labels|summary" $W/igm2_${MODE}_rd2d.log
echo "== per-Fab sub-class readout, $MODE"; $PY $W/subclass_angle_readout.py $W/route2d_$MODE/subclass_labels.txt $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt 2>&1 | grep -v Warning
echo "== same readout, $REF (5 arms, +-25 / +-30)"; $PY $W/subclass_angle_readout.py $W/route2d_$REF/subclass_labels.txt $W/sim$REF/gt_particles.txt $W/truth_$REF.txt 2>&1 | grep -v Warning
log DONE
