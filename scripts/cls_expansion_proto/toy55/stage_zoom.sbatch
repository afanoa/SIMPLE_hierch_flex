#!/bin/bash
#SBATCH -J klem_stage_zoom
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
A=/mnt/beegfs/elmlund/afan; rm -rf $A/rnasep_klem/zoom_states_f692c; /mnt/nasapps/production/python/3.12/bin/python3 $A/routes/toy55/stage_zoom_maps.py $A/rnasep_klem/f692c $A/rnasep_klem/zoom_states_f692c; du -sh $A/rnasep_klem/zoom_states_f692c; echo STAGE DONE
