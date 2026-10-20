#!/bin/bash
#SBATCH -J toy10_merge_minsize
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; GT=$W/sim/gt_particles.txt
for ms in 0 1000 2000; do
  echo "== min state size $ms"
  echo -n '[hard level 3]   '; $PY $W/expand_families.py merge $W/lat/expand_l3b_fam1 $W/lat/klem_l3b_fam1_1 $W/lat/klem_l3b_fam1_2 $W/lat/klem_l3b_fam2_1 $W/lat/klem_l3b_fam2_2 --gt $GT --minrep 0.5 --minsize $ms | grep 'MERGED' | cut -c1-330
  echo -n '[global leaves]  '; $PY $W/expand_families.py merge $W/lat/expand_l3g_fam1 $W/lat/klem_l3g_fam1_1 $W/lat/klem_l3g_fam1_2 $W/lat/klem_l3g_fam2_1 $W/lat/klem_l3g_fam2_2 --gt $GT --minrep 0.5 --bymaxr --minsize $ms | grep 'MERGED' | cut -c1-330
  echo -n '[global, own fam]'; $PY $W/expand_families.py merge $W/lat/expand_l3g_fam1 $W/lat/klem_l3g_fam1_1 $W/lat/klem_l3g_fam1_2 $W/lat/klem_l3g_fam2_1 $W/lat/klem_l3g_fam2_2 --gt $GT --minrep 0.5 --minsize $ms | grep 'MERGED' | cut -c1-330
done
