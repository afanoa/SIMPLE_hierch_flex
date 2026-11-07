#!/bin/bash
#SBATCH -J klem_umap_leaves
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
A=/mnt/beegfs/elmlund/afan; /mnt/nasapps/production/python/3.12/bin/python3 $A/routes/toy55/umap_leaves.py $A/rnasep_klem/f692c $A/routes/toy55/figs/umap_f692c_leaf_states.png 2>&1 | grep -v Warning; echo UMAP DONE
