#!/bin/bash
#SBATCH -J klem_dd
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:05:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
/mnt/nasapps/production/python/3.12/bin/python3 /mnt/beegfs/elmlund/afan/routes/toy55/truth_dist_subset.py /mnt/beegfs/elmlund/afan/routes/toy55/truthtree_55_200k.npz ${S:?}
