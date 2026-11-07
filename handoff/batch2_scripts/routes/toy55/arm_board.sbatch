#!/bin/bash
#SBATCH -J klem_arm_board
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Is the closed arm separable at all? A klem board on family 1 of the 150k dev run at its refined poses, focus = the arm's
# travel (closed + open positions), no balance, 3 x 3 lattice, host threads; then the closed-arm occupancy of every cell.
# Also the cell-level closed-arm readout of the finished f692c family-1 boards.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:-f692dev}; O=$A/rnasep_klem/arm_board_${RTAG:-f692dev}; mkdir -p $O; log(){ echo "[$(date '+%F %T') armboard] $*"; }
CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1); $PY $W/make_arm_mask.py $W/maps160/toy_consensus.mrc $CL $O/arm_mask.mrc 2>&1 | grep -v Warning
echo "== cell-level closed-arm occupancy of the f692c family-1 boards (top cells)"; for d in $A/rnasep_klem/f692c/fam1/klem_fb $A/rnasep_klem/f692c/fam1/regionmotion_1/board $A/rnasep_klem/f692c/fam1/regionmotion_3/board; do echo "-- $d"; $PY $W/closed_arm_cells.py $d --top 4 2>&1 | grep -v Warning | tail -5 | cut -c1-160; done
PS=$(ls -td $R/fam1/refine/*_refine3D_states | head -1)/p.simple; SUB=$R/expand/family1_subset.txt
for BAL in 0 0.5; do
  log "arm-focused board, balance $BAL, 3x3 lattice, family 1 ($(grep -c . $SUB) particles) at the refined poses"
  N1=3 N2=3 NFREE=0 BAL=$BAL NTHR=48 CTFABS=1 PSCALE=1 SERVER_LIST= NPARTS= FOCUSMASK=$O/arm_mask.mrc META=$R/fam1/meta_fb.npz bash $W/fboard.sh $O/board_bal$BAL $PS $SUB "12 12 10 10 8 8 8 8 8" "1.0 1.0 0.8 0.8 0.6 0.6 0.5 0.5 0.5" 2>&1 | grep -v Warning | tail -3
  echo "== closed-arm occupancy per cell (balance $BAL)"; $PY $W/closed_arm_cells.py $O/board_bal$BAL --top 9 2>&1 | grep -v Warning | tail -12 | cut -c1-160
done
log DONE
