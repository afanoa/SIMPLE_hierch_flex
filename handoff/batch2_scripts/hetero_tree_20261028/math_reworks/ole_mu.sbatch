#!/bin/bash
#SBATCH -J mrw_olemu
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
# OLE landscape scoring of the mu arm against head and bc (the ole_headbc.sbatch recipe): leaf tables, the merge-up's units, and
# landscape_score.py on the head root embedding over the leaves
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
PY=/mnt/nasapps/production/python/3.12/bin/python3; export PYTHONDONTWRITEBYTECODE=1
{ for a in head bc mu; do test -e gates/${a}_ole/1_hetero_tree/hetero_tree_leaves.txt || continue; echo "## $a ole: $(grep -c . gates/${a}_ole/1_hetero_tree/hetero_tree_leaves.txt) leaves, $(grep -h finished gates/${a}_ole/run_manifest.txt)"; awk '$2 != "junk" {s+=$3; if ($3>m) m=$3} END {printf "   largest leaf %d of %d = %.3f\n", m, s, m/s}' gates/${a}_ole/1_hetero_tree/hetero_tree_leaves.txt; grep -E 'HETERO_TREE node [a-z0-9]*: (split group|states)' gates/${a}_ole/tree.log | cut -c1-140; test -e gates/${a}_ole/1_hetero_tree/mergeup/merged_leaves.txt && echo "   merge-up units: $(grep -vc '^#' gates/${a}_ole/1_hetero_tree/mergeup/merged_leaves.txt)"; done
  echo '## landscape (leaves): head vs bc, mu; embedding = head root flex coordinates'
  $PY -I /mnt/beegfs/elmlund/afan/hetero_tree_20261028/ole_r14/tools/landscape_score.py --emb gates/head_ole/1_hetero_tree/node_r/flex/flex_pca_coordinates.txt --leaves head=gates/head_ole/1_hetero_tree bc=gates/bc_ole/1_hetero_tree mu=gates/mu_ole/1_hetero_tree 2>&1 | tail -40
} > gates/ole_mu.txt 2>&1
