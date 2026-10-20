#!/bin/bash
#SBATCH -J klem_oracle
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
W=/mnt/beegfs/elmlund/afan/routes/toy55
for d in region_1_refined_lp20 region_1_sig2_lp20 region_1_sig5_lp20 region_1; do echo "== $d"; /mnt/nasapps/production/python/3.12/bin/python3 $W/oracle_clades.py $W/armllr_55_200k/$d/arm_llr.npz $W/truthtree_55_200k.npz 2>&1 | grep -v Warning; done
