#!/bin/bash
#SBATCH -J klem_treesum
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# leaves of finished data-derived trees with their truth composition (env LIST = "MODE:TAG ...")
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
for mt in $LIST; do M=${mt%%:*}; T=${mt#*:}; echo "######## $M $T"
  for H in A B; do $PY $W/tree_leaves.py $K/synth_${M}_$H $T /tmp/leaves_${M}_${T}_$H.txt 2>&1 | grep -v Warning; $PY $W/dtree_compose.py /tmp/leaves_${M}_${T}_$H.txt $W/sim${M}_$H/gt_particles.txt 3 2>&1 | grep -v Warning | sed "s/^/   $H /"; done; done
