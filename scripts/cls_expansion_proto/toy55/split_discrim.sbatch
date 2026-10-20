#!/bin/bash
#SBATCH -J klem_splitdiscrim
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# what separates the two families of each tree (env MODES)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
for M in ${MODES:-null_50k_A null_50k_B two_50k_A two_50k_B}; do R=$A/rnasep_klem/synth_$M; echo "== $M: family 1 vs family 2"
  $PY $W/split_discrim.py $W/sim$M/all.txt $W/lat_$M/ptcl3D_r0.txt $R/expand/family1_subset.txt $R/expand/family2_subset.txt $W/sim$M/gt_particles.txt 2>&1 | grep -v Warning; done
