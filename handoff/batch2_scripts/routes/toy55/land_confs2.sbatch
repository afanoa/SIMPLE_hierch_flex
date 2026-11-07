#!/bin/bash
#SBATCH -J land_confs2
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:40:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_confs2.%j.out
# arm readouts of the deposited apo / holo conformers and of the rigid-swing maps (explicit file lists)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; O=$A/rnasep_klem/landscape/confs
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_confs.txt $T/maps160/docked_all/conf*_m160.mrc $T/maps160/holo_docked/conf*_m160.mrc $T/maps160/swing30/conf2*_swing_m160.mrc $T/maps160/swing60/conf4*_swing_m160.mrc 2>&1 | grep --line-buffered -v "Warning\|lp = lambda" | tail -2
echo "[$(date '+%F %T')] land_confs2 DONE"
