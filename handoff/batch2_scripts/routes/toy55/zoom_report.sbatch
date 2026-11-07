#!/bin/bash
#SBATCH -J zoom_report
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# arm-zoom report per half: sub-states, closed-arm readout, share of the earlier closed sub-state particles (env TAG, SELN = the zoom dir's selection name, e.g. node_r_2)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; T=$K/dtree_rnasep692_${TAG:-gs1}
for H in A B; do Z=$T/armzoom_${SELN:?}_$H; echo "== half $H"; grep -a "supported sub-states\|FAILED" $K/arm_zoom.*.out 2>/dev/null | grep "$H\]" | tail -1 | cut -c1-200
  cat $Z/arm_readout.txt 2>/dev/null | grep -v "^#" | sed "s#$Z/level1/states_dir/##" | cut -c1-120
  [ -s $Z/level1/states_clusters.txt ] && $PY $W/proxy_in_labels.py $K/f692c/fam1/zoomregionmotion_1_s05/level1/states_clusters.txt 1:3 $W/simrnasep692_$H/meta_rows.txt $Z/level1/states_clusters.txt "zoom sub-states $H"; done
