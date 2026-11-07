#!/bin/bash
#SBATCH -J ab_read
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A.out
# the state-label A/B's readout: production stage B, the argmin and the marg arms per node (f2r_ab_readout.py)
set -uo pipefail
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
G=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo PYTHONDONTWRITEBYTECODE=1
TR="closed:$G/bc_r131/fixture_rows.txt:/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt"
for k in r131:bc_r131/1_hetero_tree/node_r1:$TR:fast_r131_r1 qtoynull:bc_qtoynull/1_hetero_tree/node_r:qtoy:/mnt/beegfs/elmlund/afan/routes/toy55/simqtoynull/truth_qtoy.txt:fast_qtoynull_r olenull:bc_olenull/1_hetero_tree/node_r:none:fast_olenull_r; do
  key=${k%%:*}; rest=${k#*:}; nd=${rest%%:*}; rest=${rest#*:}; dk=${rest##*:}; truth=${rest%:*}
  N=$G/$nd; P=$L/ab_sm/${key}_prod; mkdir -p $P
  ln -sfn $L/dumps/$dk/b_ptcl3D.txt $P/b_ptcl3D.txt; ln -sfn $N/split/split_children.txt $P/split_children.txt
  echo "== $key"
  $PY $L/f2r_ab_readout.py --node $N --truth "$truth" --arm production=$P --arm argmin=$L/ab_sm/${key}_argmin --arm marg=$L/ab_sm/${key}_marg 2>&1
done
echo done
