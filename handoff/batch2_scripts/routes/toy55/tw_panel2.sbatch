#!/bin/bash
#SBATCH -J klem_treewiener2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# tree-weighted reconstruction test (tree_wiener.py) on the synthetic panel runs: leaf maps vs truth with own filter / star / tree /
# two-region tree / shuffled tree / unstructured covariance. env SETS (swing30hs prodihs two_50k null_50k)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
for t in ${SETS:-swing30hs prodihs two_50k null_50k}; do
  P=$K_/panel_$t; echo "== $t"
  $PY $W/tree_wiener2.py $K_/treewiener2/$t --rdir $P/1_refine3D_states --labels $P/refined_states.txt --tree $P/smerge/merge_tree.txt --gt $W/sim$t/gt_particles.txt --truth $W/sim$t/maps,$W/maps160,$W/maps160/docked --mskdiam 220 2>&1 | grep -v Warning
done
echo DONE
