#!/bin/bash
#SBATCH -J junk_check
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# dry run of the root-child junk flag on the roots of finished trees (MODE:TAG list in TREES, colons inside, spaces between)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
for t in ${TREES:?}; do M=${t%%:*}; T=${t#*:}; echo "== $M $T root"; $PY $W/junk_children.py $K/synth_${M}_A/$T/r $K/synth_${M}_B/$T/r r 0.5 $W/sim${M}_A/gt_particles.txt $W/sim${M}_B/gt_particles.txt 2>&1 | grep -v Warning; done
