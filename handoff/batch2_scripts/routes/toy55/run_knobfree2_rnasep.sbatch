#!/bin/bash
#SBATCH -J rnasep_knobfree2
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 16:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Overnight driver of the truth-free node recipe on family 1: the regions of the family-1 node derived from its own
# feedback board (derive_regions.py), one region_board job per region (parallel), then the landscape count over the
# region boards and the identification summary. env: MAXREG (4)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=$A/rnasep_klem/expand1; P=$E/fam1_pipeline; KD=$P/klem_fb; D=$A/rnasep_klem/knobfree; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') kf] $*"; }
log "regions of the family-1 node from its own feedback board $KD"; $PY $W/derive_regions2.py $KD $D/regions2 --ref $KD/closed_focal/focal_mask_p12.mrc --maxreg ${MAXREG:-5} 2>&1 | grep -v Warning
R=$(ls $D/regions2/region_*.mrc 2>/dev/null | wc -l); [ "$R" -ge 1 ] || { log "no regions"; exit 1; }; ids=""
for k in $(seq 1 $R); do j=$(cd $W && MASKFILE=$D/regions2/region_$k.mrc REG=p$k sbatch --export=ALL --parsable $W/region_board_rnasep.sbatch); ids="$ids $j"; log "region $k -> job $j"; done
n=0; while [ $n -lt 600 ]; do run=0; for j in $ids; do squeue -h -j $j -o %T 2>/dev/null | grep -q "RUNNING\|PENDING" && run=1; done; [ $run = 1 ] || break; sleep 120; n=$((n+1)); done
log "region boards finished: $ids"
for j in $ids; do f=$A/rnasep_klem/rnasep_region_board.$j.out; echo "===== $f"; grep -a "rb[0-9]*\]\|identified\|^ *[0-9]* *[+-][0-9.]* *[0-9]* *\(identified\|-\)\|VISIBILITY\|^  region\|^  deposited\|states ->\|state [0-9]*:" $f | cut -c1-200; done
log "landscape count over the region boards (level 1)"; $PY $W/count_states.py $(for k in $(seq 1 $R); do echo $D/region_p$k/board; done) --sim 0.6,0.7,0.8 2>&1 | grep -v Warning | tail -8
log DONE
