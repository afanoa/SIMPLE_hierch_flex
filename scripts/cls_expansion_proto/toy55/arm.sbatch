#!/bin/bash
#SBATCH -J klem_arm
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}
CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1); echo "closed reference: $CL"
echo "== consolidated 0.7 states"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $R/consolidated_0.7/arm_states.txt $R/consolidated_0.7/fam0/regionall/level1/states_dir 2>&1 | grep -v Warning
echo "== all leaf maps (level 1 and level 2)"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $R/arm_leaf.txt $(ls -d $R/fam?/region*/level1/states_dir $R/fam?/region*/level2_s*/states_dir) 2>&1 | grep -v Warning | tail -8
echo "ARM DONE"
