#!/bin/bash
#SBATCH -J klem_eocheck
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
for M in two_50k_A null_50k_A; do R=$A/rnasep_klem/synth_$M; echo "== $M"
  $PY $W/eo_check.py $W/sim$M/gt_particles.txt $R/fam1/refined_ptcl3D_print.txt,$R/fam2/refined_ptcl3D_print.txt $R/expand/family1_subset.txt,$R/expand/family2_subset.txt 2>&1 | grep -v Warning; done
