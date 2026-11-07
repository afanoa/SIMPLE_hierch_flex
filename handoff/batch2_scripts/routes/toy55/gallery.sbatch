#!/bin/bash
#SBATCH -J klem_gallery
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}; OUTN=$(echo $RTAG | tr / _); mkdir -p $W/figs
$PY $W/volume_gallery.py $R $W/maps160/toy_consensus.mrc $W/figs/gallery_$OUTN --nrow ${NROW:-12} --sort ${SORT:-mass} --lp ${LP:-10} 2>&1 | grep -v Warning; echo "GALLERY DONE $RTAG"
