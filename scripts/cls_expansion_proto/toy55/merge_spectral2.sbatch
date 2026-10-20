#!/bin/bash
#SBATCH -J toy10_merge_spectral2
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55
$PY $W/expand_families.py merge $W/lat/expand_spectral2 $W/lat/klem_expand_spectral2_fam1 $W/lat/klem_expand_spectral2_fam2 --gt $W/sim/gt_particles.txt --minrep 0.5 | sed 's/^/[spectral2] /'
