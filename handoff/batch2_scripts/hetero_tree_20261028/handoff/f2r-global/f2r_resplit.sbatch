#!/bin/bash
#SBATCH -J resplit
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# the split step of a seeding replay again under the gate trees' split environment (SIMPLE_SPLIT_CONTRAST=maps, which the
# replays lacked and which changes state_split's membership rule), into split2/, then the label readout with its children.
# LIST lines: "NAME NODE TRUTH" (NAME a directory under ab_seed/).
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=8 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo PYTHONDONTWRITEBYTECODE=1
export SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes
read -r NAME NODE TRUTH < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
R=$L/ab_seed/$NAME
test -s $R/stageB/b_ptcl3D.txt || { echo "no stage B in $R"; exit 3; }
mkdir -p $R/split2 $R/stageB2
cd $R/split2 && cp $R/stageB/b.simple s.simple
cmd=$(grep -E 'simple_exec prg=state_split' $NODE/split/hetero_step.sh | sed -e "s#^[^ ]*/simple_exec#$B/bin/simple_exec#" -e 's# > hetero_step.log.*##')
echo "$cmd" > cmd.txt
echo "== split2 $NAME $(date +%T)"; eval "$cmd" > step.log 2>&1; echo "   exit $? $(date +%T)"
grep -E 'contrast-aware|membership by the choose rule|K\* = ' step.log | cut -c1-160
ln -sfn $R/stageB/b_ptcl3D.txt $R/stageB2/b_ptcl3D.txt; ln -sfn $R/split2/split_children.txt $R/stageB2/split_children.txt
ARM=${NAME#mm2_*_}; ARM=${NAME#*_*_}
$PY $L/f2r_ab_readout.py --node $NODE --truth $TRUTH --arm production=$R/prod --arm $ARM=$R/stageB2 2>&1 | tee $R/readout/labels2.txt
echo "done $(date)"
