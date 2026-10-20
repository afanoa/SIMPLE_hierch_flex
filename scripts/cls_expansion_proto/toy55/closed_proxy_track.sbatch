#!/bin/bash
#SBATCH -J closed_track
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# track the earlier pure closed sub-states (f692c fam1 zoomregionmotion_1_s05 sub-states 1 and 3) through the real tree (env TAG, FIN)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; TAG=${TAG:-gs1}
$PY $W/closed_proxy_track.py $K/f692c/fam1/zoomregionmotion_1_s05/level1/states_clusters.txt 1:3 $K/synth_rnasep692 $TAG $W/simrnasep692 $( [ -n "${FIN:-}" ] && echo $K/dtree_rnasep692_$TAG/$FIN ) 2>&1 | grep -v Warning
