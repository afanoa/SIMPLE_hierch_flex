#!/bin/bash
#SBATCH -J toy10_merge_anneal_g
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; CH=$(ls -d $W/lat/klem_expand_anneal_g_fam* | sort -V | tr '\n' ' ')
for m in 0.0 0.3; do echo "== minrep $m"; $PY $W/expand_families.py merge $W/lat/expand_anneal_g $CH --gt $W/sim/gt_particles.txt --minrep $m | grep 'MERGED\|best-state'; done
