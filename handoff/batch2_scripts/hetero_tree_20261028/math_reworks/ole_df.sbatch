#!/bin/bash
#SBATCH -J mrw_oledf
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
PY=/mnt/nasapps/production/python/3.12/bin/python3; export PYTHONDONTWRITEBYTECODE=1
{ for a in base sub4 df; do test -e gates/${a}_ole/1_hetero_tree/hetero_tree_leaves.txt || continue; echo "## $a ole: $(grep -c . gates/${a}_ole/1_hetero_tree/hetero_tree_leaves.txt) leaves, $(grep -h finished gates/${a}_ole/run_manifest.txt)"; awk '$2 != "junk" {s+=$3; if ($3>m) m=$3} END {printf "   largest leaf %d of %d = %.3f\n", m, s, m/s}' gates/${a}_ole/1_hetero_tree/hetero_tree_leaves.txt; grep -E 'HETERO_TREE node [a-z0-9]*: (split group|states)' gates/${a}_ole/tree.log | cut -c1-140; done
  echo "## df ole split tests"; grep -h -E "draw-free\]|weakest pair at the start|locality guard: K" gates/df_ole/1_hetero_tree/node_*/split/*.log 2>/dev/null | cut -c1-200
  echo "## landscape (leaves): base vs sub4, df; embedding = base root flex coordinates"
  $PY -I /mnt/beegfs/elmlund/afan/hetero_tree_20261028/ole_r14/tools/landscape_score.py --emb gates/base_ole/1_hetero_tree/node_r/flex/flex_pca_coordinates.txt --leaves base=gates/base_ole/1_hetero_tree sub4=gates/sub4_ole/1_hetero_tree df=gates/df_ole/1_hetero_tree 2>&1 | tail -60
} > gates/ole_df.txt 2>&1
