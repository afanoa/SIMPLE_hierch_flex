#!/bin/bash
#SBATCH -J klem_pathshape
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; R=/mnt/beegfs/elmlund/afan/rnasep_klem/rflex_s30a
for nd in r22 r11 r12 r21; do $PY $T/path_shape.py $(ls -td $R/coord_$nd/leaves/*_reconstruct3D | head -1) 2>&1 | grep -v Warning; done; echo DONE
