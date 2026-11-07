#!/bin/bash
#SBATCH -J mrw_olescore
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
PY=/mnt/nasapps/production/python/3.12/bin/python3; export PYTHONDONTWRITEBYTECODE=1
{ for a in base sub; do echo "## $a ole: $(grep -c . gates/${a}_ole/1_hetero_tree/hetero_tree_leaves.txt) leaves, $(grep -h finished gates/${a}_ole/run_manifest.txt)"; awk '$2 != "junk" {s+=$3; if ($3>m) m=$3} END {printf "   largest leaf %d of %d = %.3f\n", m, s, m/s}' gates/${a}_ole/1_hetero_tree/hetero_tree_leaves.txt; grep -E 'HETERO_TREE node [a-z0-9]*: (split group|states)' gates/${a}_ole/tree.log | cut -c1-140; done
  echo '## landscape (leaves) base vs sub, embedding = base root flex coordinates'
  $PY -I /mnt/beegfs/elmlund/afan/hetero_tree_20261028/ole_r14/tools/landscape_score.py --emb gates/base_ole/1_hetero_tree/node_r/flex/flex_pca_coordinates.txt --leaves base=gates/base_ole/1_hetero_tree sub=gates/sub_ole/1_hetero_tree 2>&1 | tail -40
  echo '## 6.4 prototype on the base_ole factorial layer'
  $PY -I pole_angle_proto.py gates/base_ole/1_hetero_tree 2>&1 | tail -30
} > gates/ole_readout.txt 2>&1
