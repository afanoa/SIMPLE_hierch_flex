#!/bin/bash
#SBATCH -J mrw_sub4read
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks; export PYTHONDONTWRITEBYTECODE=1; PY=/mnt/nasapps/production/python/3.12/bin/python3
{ for t in sub4_qtoy2_s1 sub4_r131 sub4_r131_s1 r2g_qtoy2 r2g_r131; do echo "## $t ($(grep -h finished gates/$t/run_manifest.txt))"; grep -c FAILED gates/$t/tree.log; grep -E 'HETERO_TREE node [a-z0-9]*: (split group|states)' gates/$t/tree.log | cut -c1-140; grep -h 'strongest kept pole\|folded' gates/$t/tree.log | cut -c1-200 | head -4; done
  $PY -I /mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts/qtoy_leaves.py sub4_s1=gates/sub4_qtoy2_s1/1_hetero_tree r2g_s1=gates/r2g_qtoy2/1_hetero_tree 2>&1
  $PY -I /mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts/closed_nodes2.py gates/sub4_readout sub4_s2=gates/sub4_r131/1_hetero_tree sub4_s1=gates/sub4_r131_s1/1_hetero_tree r2g_s2=gates/r2g_r131/1_hetero_tree 2>&1 | grep -v -E 'Deprecation|V = np'
  sacct -u atifao -S $(date -d '-12 hours' +%Y-%m-%dT%H:%M) -X -n -o WorkDir%120,CPUTimeRAW 2>/dev/null | awk '$1 ~ "gates/(sub4|r2g)_" {split($1,p,"/"); for(i=1;i<=length(p);i++) if (p[i] ~ "^(sub4|r2g)_") ds=p[i]; c[ds]+=$2} END {for (k in c) printf "   %s %.1f core-h\n", k, c[k]/3600}'
} > gates/sub4_readout.txt 2>&1
