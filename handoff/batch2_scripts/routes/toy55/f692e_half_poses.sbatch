#!/bin/bash
#SBATCH -J f692e_half_poses
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
$PY $W/f692e_half_poses.py $K/f692e $W/simrnasep692 $W/lat_rnasep692 $K/dtree_rnasep692_f692e/poses_famrefined 2>&1 | grep -v Warning
