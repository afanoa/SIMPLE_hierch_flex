#!/bin/bash
#SBATCH -J pose_shift
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
$PY $W/pose_shift_check.py /mnt/beegfs/elmlund/afan/rnasep_coupled/basis_span_all/meta_all.npz $K/f692e $K/f692e/consolidated_12_fullw_zoom/states.txt 2>&1 | grep -v Warning
