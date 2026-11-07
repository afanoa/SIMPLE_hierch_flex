#!/bin/bash
#SBATCH -J sevprior
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# the flex prior files for refine=state_eval (seval_prior.py), one node per task. LIST lines: "KEY NODE TRUTHSPEC". The node's projects are
# dumped with the merged build under $L/dumps/KEY; the prior goes to $L/seval_prior/KEY.txt.
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
read -r KEY NODE TRUTH < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
D=$L/dumps/$KEY; mkdir -p $D $L/seval_prior $L/cache; cd $D
dump() { [ -s $3 ] || simple_exec prg=print_project_field projfile=$1 oritype=$2 mkdir=no > $3 2> $3.err; }
dump $NODE/n.simple stk n_stk.txt
dump $NODE/n.simple ptcl3D n_ptcl3D.txt
dump $NODE/stageB/b.simple ptcl3D b_ptcl3D.txt
dump $NODE/fflex/f.simple ptcl3D f_ptcl3D.txt
dump $NODE/recB0/b.simple ptcl3D r0_ptcl3D.txt
/usr/bin/time -f "wall %e s, maxrss %M kB" $PY $L/seval_prior.py --node $NODE --dump $D --cache $L/cache/sp_$KEY.npy --out $L/seval_prior/$KEY.txt \
  --nproc $SLURM_CPUS_PER_TASK > $L/seval_prior/$KEY.log 2> $L/seval_prior/${KEY}_time.txt
echo "exit $? $(tail -1 $L/seval_prior/${KEY}_time.txt)"; grep -E 'prior:|Traceback|Error' $L/seval_prior/$KEY.log | tail -3
