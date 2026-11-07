#!/bin/bash
#SBATCH -J arm_synth55
#SBATCH -p norm
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-arm readout of every leaf and zoom state map of a synthetic f692-route run (env RTAG)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; R=/mnt/beegfs/elmlund/afan/rnasep_klem/${RTAG:?}
M=$(ls $R/fam?/region*/level1/states_dir/state_*.mrc $R/fam?/region*/level2_s*/states_dir/state_*.mrc $R/fam?/zoomregion*/level1/states_dir/state_*.mrc 2>/dev/null)
$PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $W/maps160/docked_all/conf18_emd74209_m160.mrc $R/arm_all.txt $M 2>&1 | grep -v Warning | tail -3
echo "maps $(grep -vc '^#' $R/arm_all.txt); closed-arm >= 0.8: $(awk '!/^#/ && $2 >= 0.8' $R/arm_all.txt | wc -l); >= 0.6: $(awk '!/^#/ && $2 >= 0.6' $R/arm_all.txt | wc -l); top:"; grep -v '^#' $R/arm_all.txt | sort -k2 -gr | head -6 | sed "s#$R/##" | cut -c1-150
