#!/bin/bash
#SBATCH -J mrw_peek
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
export PYTHONDONTWRITEBYTECODE=1
{ grep -E "HETERO_TREE node [a-z0-9]*: (split group|states)" gates/sub_qtoy2/tree.log | cut -c1-140
  /mnt/nasapps/production/python/3.12/bin/python3 -I /mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts/qtoy_leaves.py sub=gates/sub_qtoy2/1_hetero_tree
  echo "## root site_rank outcomes (sub vs base)"; for a in base sub; do echo "--- $a"; for d in gates/${a}_qtoy2/1_hetero_tree/node_r/site gates/${a}_qtoy2/1_hetero_tree/node_r/nsite_*; do printf "%s: " $(basename $d); tr "\n" ";" < $d/site_poles.txt 2>/dev/null; echo; done; cat gates/${a}_qtoy2/1_hetero_tree/node_r/sites/chosen_site.txt 2>/dev/null; grep -h "strongest\|seed site\|seeding" gates/${a}_qtoy2/tree.log | grep "node r:" | head -3 | cut -c1-160; done
} > peek_out.txt 2>&1
