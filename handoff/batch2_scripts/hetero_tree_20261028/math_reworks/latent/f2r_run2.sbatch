#!/bin/bash
#SBATCH -J f2r2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# f2r_run.sbatch with the node's dumps made here when missing (under a lock, so array tasks sharing a KEY do not race) and the
# script chosen by SCRIPT (f2r_calib.py | f2r_fuse.py); the run goes to $L/runs/KEY/<script>_TAG. LIST lines as f2r_run.sbatch:
# "KEY NODE TRUTH TAG STAGE [args...]".
set -uo pipefail
: ${LIST:?}
SCRIPT=${SCRIPT:-f2r_calib.py}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim2/build_vd/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
read -r KEY NODE TRUTH TAG STAGE EXTRA < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST | awk '{k=$1; n=$2; t=$3; g=$4; s=$5; $1=$2=$3=$4=$5=""; print k, n, t, g, s, $0}')
test -d $NODE || { echo "no node $NODE"; exit 3; }
D=$L/dumps/$KEY; mkdir -p $D $L/cache $L/runs
(
  flock 9
  cd $D
  dump() { [ -s $3 ] || simple_exec prg=print_project_field projfile=$1 oritype=$2 mkdir=no > $3 2> $3.err; }
  dump $NODE/n.simple stk n_stk.txt
  dump $NODE/n.simple ptcl3D n_ptcl3D.txt
  dump $NODE/stageB/b.simple ptcl3D b_ptcl3D.txt
  dump $NODE/fflex/f.simple ptcl3D f_ptcl3D.txt
  dump $NODE/recB0/b.simple ptcl3D r0_ptcl3D.txt
) 9> $D/.lock
for f in n_stk n_ptcl3D b_ptcl3D f_ptcl3D r0_ptcl3D; do test -s $D/$f.txt || { echo "no dump $D/$f.txt"; exit 3; }; done
echo "dumps $KEY: n $(wc -l < $D/n_ptcl3D.txt) b $(wc -l < $D/b_ptcl3D.txt) f $(wc -l < $D/f_ptcl3D.txt) r0 $(wc -l < $D/r0_ptcl3D.txt) rows; sites $(ls $NODE/sites/site_[0-9][0-9].mrc 2>/dev/null | wc -l)"
if [ "$STAGE" = A ]; then
  DA=$L/dumps/${KEY}_stA; mkdir -p $DA; cd $DA
  flock $DA/.lock bash -c "[ -s a_ptcl3D.txt ] || simple_exec prg=print_project_field projfile=$NODE/stageA/r.simple oritype=ptcl3D mkdir=no > a_ptcl3D.txt 2> a_ptcl3D.txt.err"
  for f in n_stk n_ptcl3D f_ptcl3D r0_ptcl3D; do ln -sfn $D/$f.txt $DA/$f.txt; done
  ln -sfn $DA/a_ptcl3D.txt $DA/b_ptcl3D.txt
  D=$DA
fi
O=$L/runs/$KEY/${SCRIPT%.py}_$TAG; mkdir -p $O; cd $O
/usr/bin/time -f "wall %e s, maxrss %M kB" $PY $L/$SCRIPT --node $NODE --dump $D --out $O --cache $L/cache/f2r_${KEY}_${TAG}.npy \
  --truth $TRUTH --nproc $SLURM_CPUS_PER_TASK --cmp $L/dumps/$KEY/b_ptcl3D.txt $EXTRA > $O/f2r.log 2> $O/f2r_time.txt
echo "exit $? $(tail -1 $O/f2r_time.txt)"; tail -40 $O/f2r.log
