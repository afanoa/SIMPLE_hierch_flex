#!/bin/bash
#SBATCH -J mrw_proto
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
/mnt/nasapps/production/python/3.12/bin/python3 -I subspace_proto.py $(for d in dump_sites/*/; do test -e $d/site_gram.txt && echo $d; done) > proto_out.txt 2>&1
