#!/bin/bash
#SBATCH -J klem_truthtree
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
W=/mnt/beegfs/elmlund/afan/routes/toy55
/mnt/nasapps/production/python/3.12/bin/python3 $W/truth_dendrogram.py $W/maps160/docked_all $W/sim${MODE:-55_200k}/gt_particles.txt $W/truthtree_${MODE:-55_200k}.npz ${LP:-10} 2>&1 | grep -v Warning
