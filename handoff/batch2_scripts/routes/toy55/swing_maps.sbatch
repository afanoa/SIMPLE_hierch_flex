#!/bin/bash
#SBATCH -J klem_swingmaps
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the swing maps only (geometry check before the simulation): make_swing_maps.py into maps160/<TAG>
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; TAG=${TAG:-swing30}; M=$W/maps160/ah_norm
$PY $W/make_swing_maps.py $M/conf18_emd74209_m160.mrc $W/maps160/variants/domain_mask.mrc $W/maps160/$TAG --open $(ls $M/conf27_*_m160.mrc $M/conf28_*_m160.mrc $M/conf48_*_m160.mrc $M/conf49_*_m160.mrc | tr '\n' ',' | sed 's/,$//') --travel ${TRAVEL:-30} --nsteps ${NSTEPS:-31} 2>&1 | grep -v Warning
head -4 $W/maps160/$TAG/truth_swing.txt; tail -2 $W/maps160/$TAG/truth_swing.txt
