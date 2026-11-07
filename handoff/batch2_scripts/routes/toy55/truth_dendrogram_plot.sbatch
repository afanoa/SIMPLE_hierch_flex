#!/bin/bash
#SBATCH -J klem_truthplot
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
W=/mnt/beegfs/elmlund/afan/routes/toy55
/mnt/nasapps/production/python/3.12/bin/python3 $W/truth_dendrogram_plot.py $W/truthtree_${MODE:-55_200k}.npz $W/truthtree_${MODE:-55_200k}.png ${CUT:-0.40} 2>&1 | grep -v Warning
