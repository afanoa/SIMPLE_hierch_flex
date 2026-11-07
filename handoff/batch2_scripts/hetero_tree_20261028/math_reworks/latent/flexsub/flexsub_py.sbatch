#!/bin/bash
#SBATCH -J fsub_py
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/flexsub/logs/%x_%A_%a.out
# flexsub.py on one gate node: the within-state deformation at the stage-B poses and labels, corrected stacks under stk/<TAG>/<arm>.
# LIST lines: "TAG KEY NODE TRUTH MAXN"
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
F=$L/flexsub
PY=/mnt/nasapps/production/python/3.12/bin/python3
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONDONTWRITEBYTECODE=1
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
read -r TAG KEY NODE TRUTH MAXN < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
O=$F/runs/py_$TAG
[ -e $O/summary.json ] && { echo "EXISTS $O"; exit 3; }
mkdir -p $O $F/cache $F/stk
echo "flexsub $TAG ($KEY, maxn $MAXN) $(date)"
$PY $F/flexsub.py --node $NODE --dump $L/dumps/$KEY --out $O --stkroot $F/stk/$TAG --cache $F/cache/${TAG}_planes.npy --truth $TRUTH --maxn $MAXN 2>&1 | tee $O/flexsub.log
echo "exit ${PIPESTATUS[0]} $(date)"
