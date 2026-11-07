#!/bin/bash
#SBATCH -J klem_closed_check
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# which of the four zoomed sub-states is closed: arm readout on the NATIVE maps and on the composed ones, deposited match of both
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; Z=$A/rnasep_klem/f692e/consolidated_12_fullw_zoom; CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1); T=$A/rnasep_klem/closed_check; rm -rf "${T:?}"; mkdir -p $T/native $T/composed
for s in 01 02 03 04; do ln -sf $Z/fam0/regionall/level1/states_dir/state_$s.mrc $T/native/state_$s.mrc; ln -sf $Z/composed/state_$s.mrc $T/composed/state_$s.mrc; done
for k in native composed; do echo "== $k"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $T/arm_$k.txt $T/$k 2>&1 | grep -v Warning | grep -v "^closed-arm mask\|summary" | sed "s#$T/$k/##"; $PY $W/deposited_match.py $T/$k $W/maps160/docked $W/maps160/toy_consensus.mrc --save $T/m_$k.txt 2>&1 | grep -v Warning | grep identified; awk '!/^#/{print "  state", $1, "best conformer", $2, "corr", $3}' $T/m_$k.txt; done
echo "CLOSED CHECK DONE"
