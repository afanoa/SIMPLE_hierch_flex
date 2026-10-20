#!/bin/bash
#SBATCH -J klem_sharp_kernel
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Does the final lattice kernel dilute small states? The f692c family-1 motion_1 region board rerun with the same schedule
# but the final kernel widths sharpened (sigma 0.7 -> 0.3, then 0.15); readout = closed-arm occupancy of the cells.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692c; O=$A/rnasep_klem/sharp_kernel_f692c; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') sharp] $*"; }
PS=$(ls -td $R/fam1/refine/*_refine3D_states | head -1)/p.simple; SUB=$R/expand/family1_subset.txt; MK=$R/fam1/regions_motion/region_1.mrc
echo "== reference: the chain's own motion_1 board (final sigma 0.7)"; $PY $W/closed_arm_cells.py $R/fam1/regionmotion_1/board --top 3 2>&1 | grep -v Warning | tail -4 | cut -c1-160
for SIG in 0.3 0.15; do
  log "motion_1 board, final sigmas $SIG (10 workers)"; NTHR=8 NPARTS=10 CTFABS=1 PSCALE=1 SERVER_LIST= FOCUSMASK=$MK META=$R/fam1/meta_fb.npz bash $W/fboard.sh $O/board_sig$SIG $PS $SUB "12 12 10 10 8 8 8 8 8" "1.2 1.2 1.0 1.0 0.8 0.8 $SIG $SIG $SIG" 2>&1 | grep -v Warning | tail -2
  echo "== closed-arm occupancy per cell, final sigma $SIG"; $PY $W/closed_arm_cells.py $O/board_sig$SIG --top 5 2>&1 | grep -v Warning | tail -7 | cut -c1-160
done
log DONE
