#!/bin/bash
#SBATCH -J cell_report
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# per-cell closed-arm occupancy of an arm zoom board (both halves) + the earlier closed-sub-state share + the closed-cell subset (env TAG, SELN, OCC)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; T=$K/dtree_rnasep692_${TAG:-gs1}
for H in A B; do Z=$T/armzoom_${SELN:?}_$H; echo "== half $H"; $PY $W/closed_arm_cells.py $Z/board --top 3 2>&1 | grep -v Warning | head -4
  $PY $W/cell_report.py $Z/board $K/f692c/fam1/zoomregionmotion_1_s05/level1/states_clusters.txt 1:3 $W/simrnasep692_$H/meta_rows.txt $Z/closedcells --occ ${OCC:-0.6} 2>&1 | grep -v Warning; done
