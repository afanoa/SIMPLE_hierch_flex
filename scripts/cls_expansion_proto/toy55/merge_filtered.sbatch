#!/bin/bash
#SBATCH -J klem_merge_filt
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; R=/mnt/beegfs/elmlund/afan/rnasep_klem/f692c
echo "== map-similarity tree, leaves filtered (FSC <= 12 A, >= 1000 particles)"; $PY $W/merge_states.py $R $W/maps160/toy_consensus.mrc $R/merge_tree_filt 2>&1 | grep -v Warning | grep -E "dropped|states|similarity" | tail -12
echo "== particle-overlap tree, same filter"; $PY $W/merge_by_particles.py $R $R/merge_tree_jaccard_filt 2>&1 | grep -v Warning | grep -E "dropped|states|jaccard" | tail -8
echo MERGE FILT DONE
