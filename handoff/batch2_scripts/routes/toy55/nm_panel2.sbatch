#!/bin/bash
#SBATCH -J klem_nodemasks2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# focus masks for every node of the panel trees (node_masks.py), scored against where each node's truth maps differ. env SETS
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
for t in ${SETS:-swing30hs prodihs}; do
  P=$K_/panel_$t; echo "== $t"
  $PY $W/node_masks2.py $K_/nodemasks2/$t --rdir $P/1_refine3D_states --labels $P/refined_states.txt --tree $P/smerge/merge_tree.txt --gt $W/sim$t/gt_particles.txt --truth $W/sim$t/maps --mskdiam 220
done
echo DONE
