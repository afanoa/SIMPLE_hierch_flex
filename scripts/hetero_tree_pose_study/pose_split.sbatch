#!/bin/bash
#SBATCH -J pose_split
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/pose_study/logs/%x_%j.out
# pose-change statistics of every candidate split: the gate arms, the micelle fixture trees, the real PfCRT trees, the 692k tree
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | sed -n 's/^\(SIMPLE_[A-Z0-9_]*\)=.*/\1/p'); do unset $v; done
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/dist_sdisc/build12/install
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH OMP_NUM_THREADS=1 NPROC=16 NSHUF=200 PYTHONDONTWRITEBYTECODE=1
PY=/mnt/nasapps/production/python/3.12/bin/python3
H=/mnt/beegfs/elmlund/afan/hetero_tree_20261028; G=$H/math_reworks/gates; W=$H/fixtures/pfcrt_micelle; P=/mnt/beegfs/elmlund/afan/pfcrt_core
O=$H/pose_study; cd $O
$PY $O/pose_split.py --selftest || exit 1
SET=${SET:-gates}
case $SET in
  gates)
    ARGS=""
    for arm in base sub df; do
      for ds in qtoy2 r131 r131_s1 qtoynull micnull; do ARGS="$ARGS ${arm}_${ds}=$G/${arm}_${ds}"; done
      for ds in ole olenull; do ARGS="$ARGS ${arm}_${ds}=$G/${arm}_${ds}:c2"; done
    done
    for t in mic3_order mic3_nomask mic3_umask micnull_order micnull_nomask micnull_umask mic3_fp_order mic3_gate_order mic3_mixB_order; do
      test -d $W/trees/$t && ARGS="$ARGS fx_$t=$W/trees/$t"
    done
    for t in tree_umask_s1 tree_umask2_s1; do test -d $P/$t && ARGS="$ARGS pfcrt_$t=$P/$t"; done
    $PY $O/pose_split.py $O/pose_gates $ARGS ;;
  rnasep)
    $PY $O/pose_split.py $O/pose_rnasep rnasep692=$H/integration/f692_r17d_s1 ;;
esac
