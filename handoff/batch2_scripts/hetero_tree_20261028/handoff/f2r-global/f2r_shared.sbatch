#!/bin/bash
#SBATCH -J f2rshared
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# seam-free global coordinate of one tree per array task (f2r_shared.py): the root project's dumps, every leaf node's stage-B
# dump, the final poses assembled (f2r_leafposes.py), then three embeddings: the tree-wide basis at the final poses, the root's
# modes at the final poses (control), the tree-wide basis at the input poses. LIST lines: "KEY TREE TRUTH CACHE [args]".
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo MEST_LATENT=$L MEST_DIR=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest
read -r KEY TREE TRUTH CACHE EXTRA < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST | awk '{k=$1; t=$2; r=$3; c=$4; $1=$2=$3=$4=""; print k, t, r, c, $0}')
test -d $TREE/node_r || { echo "no tree $TREE"; exit 3; }
D=$L/dumps/root_$KEY; mkdir -p $D/input $D/final
(
  flock 9
  cd $D
  dump() { [ -s $3 ] || simple_exec prg=print_project_field projfile=$1 oritype=$2 mkdir=no > $3 2> $3.err; }
  dump $TREE/node_r/n.simple stk n_stk.txt
  dump $TREE/node_r/n.simple ptcl3D n_ptcl3D.txt
  for x in $(awk '{print $1}' $TREE/hetero_tree_leaves.txt | sed 's/_.*//' | sort -u); do
    [ -e $TREE/node_$x/stageB/b.simple ] && dump $TREE/node_$x/stageB/b.simple ptcl3D leaf_${x}_ptcl3D.txt
  done
  [ -s final_ptcl3D.txt ] || $PY $L/f2r_leafposes.py $TREE n_ptcl3D.txt $D final_ptcl3D.txt
) 9> $D/.lock
echo "dumps: root $(wc -l < $D/n_ptcl3D.txt) rows, $(ls $D/leaf_*_ptcl3D.txt | wc -l) leaf dumps, final $(wc -l < $D/final_ptcl3D.txt) rows"
ln -sfn $D/n_stk.txt $D/input/n_stk.txt; ln -sfn $D/n_ptcl3D.txt $D/input/n_ptcl3D.txt
ln -sfn $D/n_stk.txt $D/final/n_stk.txt; ln -sfn $D/n_ptcl3D.txt $D/final/n_ptcl3D.txt; ln -sfn $D/final_ptcl3D.txt $D/final/b_ptcl3D.txt
for run in "all final" "root final" "all input"; do
  set -- $run; O=$L/shared/${KEY}_$1_$2; mkdir -p $O
  echo "==== $KEY basis $1 poses $2 $(date +%T)"
  $PY $L/f2r_shared.py embed --tree $TREE --node $TREE/node_r --dump $D/$2 --cache $CACHE --out $O --truth "$TRUTH" --basis $1 --poses $2 --nproc $SLURM_CPUS_PER_TASK $EXTRA > $O/run.log 2>&1
  echo "   exit $? $(date +%T)"; grep -E 'Error|Traceback' $O/run.log | head -3
  grep -E '^basis|^embedding|== \(a\)|named pair|AUROC\(d' $O/report.txt | sed 's/^\[[0-9:]*\] *//' | cut -c1-200
done
echo "done $(date)"
