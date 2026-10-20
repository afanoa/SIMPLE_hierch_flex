#!/bin/bash
#SBATCH -J klem_dtplot
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# render a data-derived tree (env MODE, TAG, TITLE)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:?}
$PY $W/dtree_plot.py $K/synth_${MODE}_A $K/synth_${MODE}_B $TAG $W/sim${MODE}_A/gt_particles.txt $W/sim${MODE}_B/gt_particles.txt $K/dtree_${MODE}_$TAG/tree.png --truth $W/truthtree_55_200k.npz --title "${TITLE:-}" 2>&1 | grep -v Warning
