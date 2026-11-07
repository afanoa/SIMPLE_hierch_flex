#!/bin/bash
#SBATCH -J klem_pooled_analysis
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# readouts of a consolidated (pooled) state set: env RTAG (run), CDIR (consolidated_<cut><tag>): gallery with FSC/rep per
# state, arm readout (closed/open occupancy), strict recovery of the 55 deposited apo conformers
set -uo pipefail; unset LEAFGLOB
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}; C=$R/${CDIR:?}; SD=$C/fam0/regionall/level1/states_dir; OUTN=$(echo ${RTAG}_${CDIR} | tr / _)
echo "== $CDIR: $(ls $SD/state_*.mrc | wc -l) pooled states"; head -40 $C/groups.txt | cut -c1-120
echo "== gallery (fsc_A per state)"; $PY $W/volume_gallery.py $C $W/maps160/toy_consensus.mrc $W/figs/gallery_$OUTN --nrow ${NROW:-8} --sort ${SORT:-mass} --lp ${LP:-7} 2>&1 | grep -v Warning | tail -2
awk '!/^#/{print $1, $2, $3, $5, $6}' $W/figs/gallery_${OUTN}_table.txt | head -60
CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1); echo "== arm readout (closed/open occupancy per state)"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $C/arm_states.txt $SD 2>&1 | grep -v Warning | grep -v mask | sed "s#$SD/##" | sort -k3 -rn | head -12; echo "closed-called: $(grep -c 'closed$' $C/arm_states.txt)"
echo "== strict recovery of 55"; S=$C/readout55; rm -rf $S; mkdir -p $S/states; for m in $SD/state_*.mrc; do ln -sf $m $S/states/$(basename $m); done; ls $S/states > $S/map_order.txt
$PY $W/deposited_match.py $S/states $W/maps160/docked_all $W/maps160/toy_consensus.mrc --save $S/matches.txt 2>&1 | grep -v Warning | tail -2; $PY $W/recover55.py $S/matches.txt --names $S/map_order.txt 2>&1 | grep -v Warning | head -3
echo "POOLED DONE $CDIR"
