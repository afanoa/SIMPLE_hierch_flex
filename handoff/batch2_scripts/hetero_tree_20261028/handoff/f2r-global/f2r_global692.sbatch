#!/bin/bash
#SBATCH -J f2rg692
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A.out
# global chart coordinate tests on the 692k trees (m33 and its sister recipes), then the cross-recipe comparisons.
set -uo pipefail
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
S=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo MEST_LATENT=$L MEST_DIR=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest
O=$L/global; mkdir -p $O
TR="closed:/nonexistent:/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt"
for k in full692_m33_s1 full692_base92_s1 full692_base80_s1; do
  T=$S/$k/1_hetero_tree
  echo "==== compose $k $(date +%T)"
  pair=""; grep -q '^r1111221_s1 ' $T/hetero_tree_leaves.txt && grep -q '^r1111221_s2 ' $T/hetero_tree_leaves.txt && pair="--pair r1111221_s1,r1111221_s2"
  $PY $L/f2r_global.py compose $T $O/$k --truth "$TR" --pairs siblings --merged $T/hetero_tree_merged.txt $pair --dsub 60000 --maxleaves 15 > $O/$k.log 2>&1; echo "   exit $? $(date +%T)"; grep -E 'Error|Traceback' $O/$k.log | head -3
done
echo "==== compare"
for pair in "full692_m33_s1 full692_base92_s1" "full692_m33_s1 full692_base80_s1" "full692_base92_s1 full692_base80_s1"; do
  set -- $pair
  $PY $L/f2r_global.py compare $O/$1 $O/$2 $1 $2 2>&1 | tee $O/compare_$1_vs_$2.txt
done
echo "done $(date)"
