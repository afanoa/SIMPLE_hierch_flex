#!/bin/bash
#SBATCH -J mest_charts
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest/logs/%x_%A_%a.out
# chart coordination and the density-gap split criterion of one finished tree per array task. LIST lines: "KEY TREE TRUTH [args]"
set -uo pipefail
: ${LIST:?}
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo MEST_LATENT=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
read -r KEY TREE TRUTH EXTRA < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST | awk '{k=$1; n=$2; t=$3; $1=$2=$3=""; print k, n, t, $0}')
test -d $TREE || { echo "no tree $TREE"; exit 3; }
O=$M/runs/charts_$KEY; mkdir -p $O; cd $O
/usr/bin/time -f "wall %e s, maxrss %M kB" $PY -I $M/mest_charts.py $TREE $O --truth $TRUTH --nproc $SLURM_CPUS_PER_TASK $EXTRA > run.log 2> time.txt
echo "exit $? $(tail -1 time.txt)"; grep -v "^\[" run.log | head -5; cat charts.txt
