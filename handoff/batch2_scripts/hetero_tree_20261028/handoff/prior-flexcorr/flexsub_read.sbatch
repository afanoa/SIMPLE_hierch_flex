#!/bin/bash
#SBATCH -J fsub_read
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/flexsub/logs/%x_%j.out
# the flexsub replays' state references (even/odd FSC per state) beside production's, and the positive-control toy's arm cut.
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
F=$L/flexsub
G=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo OMP_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
for key in ${KEYS:-qtoyflex qtoy2 qtoynull}; do
  N=$G/mm2_$key/1_hetero_tree/node_r
  args="production=$N/stageB"
  for arm in orig origs2 full perp; do d=$F/runs/mm2_${key}_$arm/stageB; [ -d $d ] && args="$args $arm=$d"; done
  echo "== $key state references"; $PY $F/flexsub_fsc.py 225.456 $args
done
N=$G/mm2_qtoyflex/1_hetero_tree/node_r
args="production=$F/runs/mm2_qtoyflex_orig/prod"
for arm in orig origs2 full perp; do d=$F/runs/mm2_qtoyflex_$arm/stageB; [ -s $d/b_ptcl3D.txt ] && args="$args $arm=$d"; done
echo "== qtoyflex arm cut"; $PY $F/flexsub_cut.py $N /mnt/beegfs/elmlund/afan/routes/toy55/simqtoyflex/truth_qtoy.txt $args
