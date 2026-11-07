#!/bin/bash
#SBATCH -J latent
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# fixed-pose latent-space state test, one node per array task. LIST lines: "KEY NODE TRUTH [latent_states.py args...]".
# The node's projects are dumped once (print_project_field) and its planes cached once under $L/cache/KEY.npy. SELFTEST=1 runs the
# synthetic check of the latent EM first.
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim2/build_vd/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
read -r KEY NODE TRUTH EXTRA < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST | awk '{k=$1; n=$2; t=$3; $1=$2=$3=""; print k, n, t, $0}')
test -d $NODE || { echo "no node $NODE"; exit 3; }
if [ "${SELFTEST:-0}" = 1 ]; then $PY $L/latent_states_v2.py --selftest --nproc $SLURM_CPUS_PER_TASK || exit 4; fi
D=$L/dumps/$KEY; mkdir -p $D $L/cache $L/runs; cd $D
dump() { [ -s $3 ] || simple_exec prg=print_project_field projfile=$1 oritype=$2 mkdir=no > $3 2> $3.err; }
dump $NODE/n.simple stk n_stk.txt
dump $NODE/n.simple ptcl3D n_ptcl3D.txt
dump $NODE/stageB/b.simple ptcl3D b_ptcl3D.txt
dump $NODE/fflex/f.simple ptcl3D f_ptcl3D.txt
dump $NODE/recB0/b.simple ptcl3D r0_ptcl3D.txt
echo "dumps $KEY: n $(wc -l < n_ptcl3D.txt) b $(wc -l < b_ptcl3D.txt) f $(wc -l < f_ptcl3D.txt) r0 $(wc -l < r0_ptcl3D.txt) rows, $(wc -l < n_stk.txt) stacks"
O=$L/runs/$KEY; mkdir -p $O; cd $O
/usr/bin/time -f "wall %e s, maxrss %M kB" $PY $L/latent_states_v2.py --node $NODE --dump $D --out $O --cache $L/cache/$KEY.npy \
  --truth $TRUTH --nproc $SLURM_CPUS_PER_TASK $EXTRA > run.log 2> time.txt
echo "exit $? $(tail -1 time.txt)"; tail -45 run.log
