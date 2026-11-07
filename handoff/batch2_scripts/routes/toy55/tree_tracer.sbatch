#!/bin/bash
#SBATCH -J klem_tracer
#SBATCH -p norm
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
for rt in ${RTS:-f692rf}; do echo "== rflex_$rt"; $PY $T/tree_tracer.py $K_/rflex_$rt $K_/f692r/coordtree_fam1/units.txt 16,15 2>&1 | grep -v Warning; done; echo DONE
