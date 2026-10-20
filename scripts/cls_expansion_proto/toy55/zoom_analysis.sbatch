#!/bin/bash
#SBATCH -J klem_zoom_analysis
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# after the zoom stage of a run (env RTAG): merge trees, gallery, arm readout and strict 55-conformer recovery over the ZOOMED leaves
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}; export LEAFGLOB='zoomregion*'
CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1)
echo "== zoomed leaves: $(ls $R/fam?/zoomregion*/level1/states_dir/state_*.mrc | wc -l) sub-state maps from $(ls -d $R/fam?/zoomregion* | wc -l) zooms"
echo "== map-similarity tree (FSC-filtered)"; $PY $W/merge_states.py $R $W/maps160/toy_consensus.mrc $R/merge_tree_zoom 2>&1 | grep -v Warning | grep -E "similarity|states" | tail -7
echo "== particle-overlap tree (FSC-filtered)"; $PY $W/merge_by_particles.py $R $R/merge_tree_zoom_jaccard 2>&1 | grep -v Warning | grep -E "jaccard|states" | tail -6
echo "== arm readout, top 8 by closed occupancy"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $R/arm_zoom.txt $(ls -d $R/fam?/zoomregion*/level1/states_dir) 2>&1 | grep -v Warning | grep -v mask | sed "s#$R/##" | sort -k3 -rn | head -8
grep -c "closed$" $R/arm_zoom.txt | sed 's/^/closed-called maps: /'
echo "== strict recovery of 55 (zoomed leaves)"; S=$R/readout55_zoom; rm -rf $S; mkdir -p $S/states
for d in $R/fam?/zoomregion*/level1/states_dir; do t=$(echo $d | sed "s#$R/##; s#/level1/states_dir##; s#/#_#g"); for m in $d/state_*.mrc; do ln -sf $m $S/states/state_${t}_$(basename $m .mrc | sed "s/^state_//").mrc; done; done
ls $S/states > $S/map_order.txt; $PY $W/deposited_match.py $S/states $W/maps160/docked_all $W/maps160/toy_consensus.mrc --save $S/matches.txt > /dev/null 2>&1; $PY $W/recover55.py $S/matches.txt --names $S/map_order.txt 2>&1 | grep -v Warning | head -1
echo "== gallery"; $PY $W/volume_gallery.py $R $W/maps160/toy_consensus.mrc $W/figs/gallery_${RTAG}_zoom --nrow 12 --lp 7 --sort rep 2>&1 | grep -v Warning | tail -1
echo "ZOOM ANALYSIS DONE"
