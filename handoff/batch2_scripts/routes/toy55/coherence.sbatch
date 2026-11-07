#!/bin/bash
#SBATCH -J klem_coherence
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; RTAG=${RTAG:?}; MODE=${MODE:?}
for C in ${CDIRS:-consolidated_12_rule_12 consolidated_12_rule_auto}; do echo "== $RTAG/$C"; $PY $W/state_coherence.py $A/rnasep_klem/$RTAG/$C $W/sim$MODE/gt_particles.txt $W/maps160/docked_all $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning; done; echo "COHERENCE DONE"
