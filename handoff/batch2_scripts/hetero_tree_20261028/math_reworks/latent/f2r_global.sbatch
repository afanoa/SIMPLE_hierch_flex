#!/bin/bash
#SBATCH -J f2rglob
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A.out
# global chart coordinate tests (f2r_global.py): compose on the gate trees, then the cross-tree comparisons; one job in sequence.
set -uo pipefail
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
G=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo MEST_LATENT=$L MEST_DIR=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest
O=$L/global; mkdir -p $O
TR="closed:$G/mm2_r131/fixture_rows.txt:/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt"
TQ=qtoy:/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/truth_qtoy.txt
TF=qtoy:/mnt/beegfs/elmlund/afan/routes/toy55/simqtoyflex/truth_qtoy.txt
for spec in "mm2_r131 $TR" "mm2_r131_s1 $TR" "pb15_r131 $TR" "mm2_qtoy2 $TQ" "pb15_qtoy2 $TQ" "mm2_qtoyflex $TF"; do
  set -- $spec; k=$1; t=$2
  echo "==== compose $k $(date +%T)"
  $PY $L/f2r_global.py compose $G/$k/1_hetero_tree $O/$k --truth "$t" > $O/$k.log 2>&1; echo "   exit $?"; grep -vE '^\[' $O/$k.log | head -3
done
echo "==== compare"
for pair in "mm2_r131 mm2_r131_s1" "mm2_r131 pb15_r131" "mm2_r131_s1 pb15_r131" "mm2_qtoy2 pb15_qtoy2"; do
  set -- $pair
  $PY $L/f2r_global.py compare $O/$1 $O/$2 $1 $2 2>&1 | tee $O/compare_$1_vs_$2.txt
done
echo "done $(date)"
