#!/bin/bash
#SBATCH -J klem_soft_check
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; d=$A/rnasep_klem/f692c/fam1/regionmotion_1
/mnt/nasapps/production/python/3.12/bin/python3 $W/soft_state_maps.py $d/board $d/level1/clusters.txt $W/maps160/toy_consensus.mrc /tmp/soft_check_$$ 2>&1 | grep -v Warning | head -3; echo CHECK DONE
