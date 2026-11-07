#!/bin/bash
#SBATCH -J mest_prior
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest/logs/%x_%j.out
# final-map estimators (tree-prior variants of tree_wiener) on a finished tree's raw leaf halves; env: TREE (1_hetero_tree dir),
# OUT; optional TRUTH/CSV/MAPS (qtoy scoring via mest_score.py), CLOSED=1 (RNase P closed-arm readout)
set -uo pipefail
: ${TREE:?} ${OUT:?}
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
ARGS=""
[ -n "${TRUTH:-}" ] && ARGS="$ARGS --truth $TRUTH --csv $CSV --maps $MAPS"
[ "${CLOSED:-0}" = 1 ] && ARGS="$ARGS --closed"
mkdir -p $OUT
/usr/bin/time -f "wall %e s, maxrss %M kB" $PY -I $M/mest_prior.py $TREE $OUT $ARGS --ests emp2,topo2,depth2,anc2 --closed_top 8 --halves_only > $OUT/run.log 2> $OUT/time.txt
echo "exit $? $(tail -1 $OUT/time.txt)"
cat $OUT/run.log
test -s $OUT/score/mest_score.txt && cat $OUT/score/mest_score.txt
