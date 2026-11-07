#!/bin/bash
#SBATCH -J rnasep_score_maps
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; P=/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline
echo '#### ten matched-state maps (labelled by conformer)'; $PY $W/score_maps.py $P/klem_fb/matched_rec/maps --labels 1,10,18,27,28,37,46,48,49,51
echo '#### the ten refined states of family 1 (unlabelled)'; $PY $W/score_maps.py $P/maps
