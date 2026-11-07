#!/bin/bash
#SBATCH -J ab_sm
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# stage-B A/B of the state label rule on a gate node: the node's own stage-B refine3D command (stageB/hetero_step.sh) rerun from
# the stage-B input project (recB/b.simple: stage-A poses, folded labels) with the merge_sm build, ARM argmin (as the tree) or
# marg (prob_state_marg=yes), then the node's split step (split/hetero_step.sh) on the result and the ptcl3D dump for the
# readout. The master submits its part jobs itself as the tree's did. LIST lines: "KEY NODE ARM".
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_sm/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=8 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
read -r KEY NODE ARM < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
test -d $NODE/recB || { echo "no recB in $NODE"; exit 3; }
O=$L/ab_sm/${KEY}_$ARM
[ -e $O ] && { echo "EXISTS $O"; exit 3; }
mkdir -p $O && cd $O
cp $NODE/recB/b.simple b.simple
cmd=$(grep 'prg=refine3D' $NODE/stageB/hetero_step.sh | sed -e "s#^[^ ]*/simple_exec#$B/bin/simple_exec#" -e 's# > hetero_step.log.*##')
[ "$ARM" = marg ] && cmd="$cmd prob_state_marg=yes"
echo "$cmd" > stageB_cmd.txt
echo "stage B ($ARM): $(date)"
eval "$cmd" > stageB.log 2>&1; echo "refine3D exit $?"
grep -c "NORMAL STOP\|>>> DONE\|CONVERGED" stageB.log | sed 's/^/stage B done lines: /'
cp b.simple s.simple
scmd=$(grep 'prg=state_split' $NODE/split/hetero_step.sh | sed -e "s#^[^ ]*/simple_exec#$B/bin/simple_exec#" -e 's# > hetero_step.log.*##')
echo "$scmd" > split_cmd.txt
eval "$scmd" > split.log 2>&1; echo "state_split exit $?"
echo "split children: $(sort split_children.txt 2>/dev/null | uniq -c | tr '\n' ' ')"
simple_exec prg=print_project_field projfile=b.simple oritype=ptcl3D mkdir=no > b_ptcl3D.txt 2> b_ptcl3D.err
echo "dump rows $(wc -l < b_ptcl3D.txt); $(date)"
