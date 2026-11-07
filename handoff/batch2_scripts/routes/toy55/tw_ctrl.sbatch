#!/bin/bash
#SBATCH -J klem_treewiener_ctrl
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# tree_wiener2.py on the two control sets with the truth maps taken from maps160/ah_norm, then maps160/docked (frame check decides)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
for t in two_50k null_50k; do for d in ah_norm docked; do
  P=$K_/panel_$t; echo "== $t, truth $d"
  $PY $W/tree_wiener2.py $K_/treewiener2/${t}_$d --rdir $P/1_refine3D_states --labels $P/refined_states.txt --tree $P/smerge/merge_tree.txt --gt $W/sim$t/gt_particles.txt --truth $W/maps160/$d --mskdiam 220
done; done
echo DONE
