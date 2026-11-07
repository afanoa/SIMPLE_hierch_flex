#!/bin/bash
#SBATCH -J mest_estep
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest/logs/%x_%A_%a.out
# label-level M-estimation formulations at frozen stage-B poses, one node per array task. LIST lines: "KEY NODE TRUTH [args]";
# the node's dumps and plane cache are the latent session's (read only; a missing cache is written under $M/cache).
set -uo pipefail
: ${LIST:?}
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo MEST_LATENT=$L
read -r KEY NODE TRUTH EXTRA < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST | awk '{k=$1; n=$2; t=$3; $1=$2=$3=""; print k, n, t, $0}')
test -d $NODE || { echo "no node $NODE"; exit 3; }
D=$L/dumps/$KEY; test -s $D/b_ptcl3D.txt || { echo "no dumps at $D"; exit 3; }
C=$L/cache/$KEY.npy; test -s $C || { mkdir -p $M/cache; C=$M/cache/$KEY.npy; }
O=$M/runs/estep_$KEY; mkdir -p $O; cd $O
/usr/bin/time -f "wall %e s, maxrss %M kB" $PY -I $M/mest_estep.py --node $NODE --dump $D --out $O --cache $C --truth $TRUTH --nproc $SLURM_CPUS_PER_TASK $EXTRA > run.log 2> time.txt
echo "exit $? $(tail -1 time.txt)"; tail -40 run.log
