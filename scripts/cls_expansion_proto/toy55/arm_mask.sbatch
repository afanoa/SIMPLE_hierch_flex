#!/bin/bash
#SBATCH -J arm_mask
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# arm-travel mask for the real tree's arm-focused zooms (env TAG)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; T=$K/dtree_rnasep692_${TAG:-gs1}
$PY $W/make_arm_mask.py $A/rnasep_states/25_refine3D_states/rec_final_state01.mrc $W/maps160/docked_all/conf18_emd74209_m160.mrc $T/arm_mask.mrc 2>&1 | grep -v Warning
