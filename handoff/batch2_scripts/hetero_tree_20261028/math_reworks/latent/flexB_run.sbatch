#!/bin/bash
#SBATCH -J flexB
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# flex_pca rerun at the node's stage-B poses (frozen): the tree's own fflex call (fflex/hetero_step.sh, same binary, mean and mask)
# on a copy of stageB/b.simple, then the latent test of its modes (flex_refit3_test.py, orig arm, without and with site factors).
# LIST lines: "KEY NODE TRUTH [args...]"; outputs in $L/flexB/KEY.
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
BC=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/build_bc/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$BC/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$BC
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
read -r KEY NODE TRUTH EXTRA < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST | awk '{k=$1; n=$2; t=$3; $1=$2=$3=""; print k, n, t, $0}')
D=$L/flexB/$KEY; mkdir -p $D/stderrout; cd $D
if [ ! -s FLEXB_DONE ]; then
  cp $NODE/stageB/b.simple f.simple
  CMD=$(grep 'prg=flex_pca' $NODE/fflex/hetero_step.sh | sed 's/ *> *hetero_step.log 2>&1 *$//')
  echo "$CMD" > flexB_cmd.txt
  OMP_NUM_THREADS=8 /usr/bin/time -f "flex_pca wall %e s, maxrss %M kB" $CMD > flexB.log 2> flexB_time.txt
  rc=$?; echo "flex_pca exit $rc: $(tail -1 flexB_time.txt), $(ls flex_pca_pc[0-9][0-9][0-9].mrc 2>/dev/null | wc -l) modes"
  [ $rc = 0 ] || exit $rc
  echo ok > FLEXB_DONE
fi
for sf in no yes; do
  OMP_NUM_THREADS=1 /usr/bin/time -f "eval wall %e s, maxrss %M kB" $PY $L/flex_refit3_test.py --node $NODE --dump $L/dumps/$KEY --out $D/eval_sf$sf \
    --cache $L/cache/fp_$KEY.npy --truth $TRUTH --nproc $SLURM_CPUS_PER_TASK $EXTRA --arm orig --sf $sf --pcdir $D > eval_sf$sf.log 2> eval_sf${sf}_time.txt
  echo "eval sf $sf exit $? $(tail -1 eval_sf${sf}_time.txt)"
done
