#!/bin/bash
#SBATCH -J fsub_rep
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/flexsub/logs/%x_%A_%a.out
# a gate node's stage B and split replayed on a relocated project: the node's recB input project (stage A's poses and folded
# labels) remapped to the arm's stack root by update_project, then the node's saved recB, stage B and split calls (merge_master2
# build, the tree's seed and split environment). Readouts against the truth and production (f2r_ab_readout.py, flexsub_cut.py).
# LIST lines: "KEY NODE ARM SEED OLDROOT NEWROOT TRUTH KINDTRUTH"
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
F=$L/flexsub
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=8 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo PYTHONDONTWRITEBYTECODE=1
export SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes
read -r KEY NODE ARM SEED OLDROOT NEWROOT TRUTH KINDTRUTH < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
export SIMPLE_SEED=$SEED
R=$F/runs/${KEY}_${ARM}
[ -e $R ] && { echo "EXISTS $R"; exit 3; }
mkdir -p $R/recB $R/stageB $R/split $R/prod $R/readout
echo "replay $KEY $ARM seed $SEED: $OLDROOT -> $NEWROOT; $(date)"
step() { # DIR SAVED SEDX: the saved step's simple_exec call with the binary and SEDX swaps, run in $R/DIR
  local d=$1 src=$2 sedx=$3 cmd
  cd $R/$d
  cmd=$(grep -E 'simple_exec prg=' $NODE/$src/hetero_step.sh | sed -e "s#^[^ ]*/simple_exec#$B/bin/simple_exec#" -e 's# > hetero_step.log.*##' -e "$sedx")
  echo "$cmd" > cmd.txt
  echo "== $d $(date +%T)"; eval "$cmd" > step.log 2>&1; echo "   exit $? $(date +%T)"
}
cp $NODE/recB/b.simple $R/recB/b0.simple
( cd $R/recB && $B/bin/simple_exec prg=update_project projfile=b0.simple ptcl_old_root=$OLDROOT ptcl_new_root=$NEWROOT \
    projfile_out=$R/recB/b.simple mkdir=no > update.log 2>&1 ); echo "update_project exit $?"
grep -E 'REMAPPED|WROTE|ERROR|error' $R/recB/update.log
test -s $R/recB/b.simple || { echo "no remapped project"; exit 4; }
simple_exec prg=print_project_field projfile=$R/recB/b.simple oritype=stk mkdir=no 2> /dev/null | grep -o 'stk=[^ ]*' | sort | uniq -c
step recB recB "s#x#x#"
cp $R/recB/b.simple $R/stageB/b.simple
step stageB stageB "s#$NODE/recB/#$R/recB/#g"
cp $R/stageB/b.simple $R/split/s.simple
step split split "s#x#x#"
grep -E 'contrast-aware|membership by the choose rule|K\* = ' $R/split/step.log | cut -c1-160
simple_exec prg=print_project_field projfile=$R/stageB/b.simple oritype=ptcl3D mkdir=no > $R/stageB/b_ptcl3D.txt 2> /dev/null
ln -sfn $R/split/split_children.txt $R/stageB/split_children.txt
ln -sfn $L/dumps/$KEY/b_ptcl3D.txt $R/prod/b_ptcl3D.txt; ln -sfn $NODE/split/split_children.txt $R/prod/split_children.txt
echo "== readout"
$PY $L/f2r_ab_readout.py --node $NODE --truth $TRUTH --arm production=$R/prod --arm $ARM=$R/stageB 2>&1 | tee $R/readout/labels.txt
if [ "$KINDTRUTH" != "none" ]; then
  echo "-- kind-collapsed truth"
  $PY $L/f2r_ab_readout.py --node $NODE --truth $KINDTRUTH --arm production=$R/prod --arm $ARM=$R/stageB 2>&1 | tee $R/readout/labels_kind.txt
  $PY $F/flexsub_cut.py $NODE ${TRUTH#qtoy:} production=$R/prod $ARM=$R/stageB 2>&1 | tee $R/readout/cut.txt
fi
echo "done $(date)"
