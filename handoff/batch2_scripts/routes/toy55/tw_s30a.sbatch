#!/bin/bash
#SBATCH -J klem_treewiener_s30a
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# tree-weighted leaf maps of rflex_s30a against truth (box-88 leaves, truth Fourier-cropped)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; R=/mnt/beegfs/elmlund/afan/rnasep_klem/rflex_s30a; P=$R/post
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
$PY $T/tree_wiener5.py $P/tw --list $P/leaves_int.txt --labels $R/final_labels.txt --tree $P/tree/merge_tree.txt --gt $T/simswing30/gt_particles.txt --truth $T/simswing30/maps --mskdiam 220 2>&1 | grep -v Warning
echo DONE
