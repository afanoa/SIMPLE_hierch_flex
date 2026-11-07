#!/bin/bash
#SBATCH -J klem_sd2
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
for H in A B; do echo "== half $H: child 1 vs child 3, conformer 49 only"
/mnt/nasapps/production/python/3.12/bin/python3 $W/split_discrim.py $W/simtwo_50k_$H/all.txt $W/lat_two_50k_$H/ptcl3D_r0.txt $K/synth_two_50k_$H/dtree5/r/children.txt:1 $K/synth_two_50k_$H/dtree5/r/children.txt:3 $W/simtwo_50k_$H/gt_particles.txt 49 2>&1 | grep -v Warning | head -12; done
