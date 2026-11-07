#!/bin/bash
#SBATCH -J f2r
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# flex -> refinement replication (f2r_calib.py), one run per array task. LIST lines: "KEY NODE TRUTH TAG STAGE [f2r_calib.py args...]".
# KEY names the node's dumps ($L/dumps/KEY, made by calib2_run.sbatch); TAG names the run ($L/runs/KEY/f2r_TAG) and its own plane
# cache ($L/cache/f2r_KEY_TAG.npy: prep_crop takes a cache whose shape matches, so another subsample must not share one).
# STAGE B runs on stage B's labels and poses; STAGE A on stage A's (stageA/r.simple dumped once into $L/dumps/KEY_stA, linked there
# as b_ptcl3D.txt beside the node's other dumps).
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim2/build_vd/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
read -r KEY NODE TRUTH TAG STAGE EXTRA < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST | awk '{k=$1; n=$2; t=$3; g=$4; s=$5; $1=$2=$3=$4=$5=""; print k, n, t, g, s, $0}')
test -d $NODE || { echo "no node $NODE"; exit 3; }
D=$L/dumps/$KEY
for f in n_stk n_ptcl3D b_ptcl3D f_ptcl3D r0_ptcl3D; do test -s $D/$f.txt || { echo "no dump $D/$f.txt"; exit 3; }; done
if [ "$STAGE" = A ]; then
  DA=$L/dumps/${KEY}_stA; mkdir -p $DA; cd $DA
  flock $DA/.lock bash -c "[ -s a_ptcl3D.txt ] || simple_exec prg=print_project_field projfile=$NODE/stageA/r.simple oritype=ptcl3D mkdir=no > a_ptcl3D.txt 2> a_ptcl3D.txt.err"
  for f in n_stk n_ptcl3D f_ptcl3D r0_ptcl3D; do ln -sfn $D/$f.txt $DA/$f.txt; done
  ln -sfn $DA/a_ptcl3D.txt $DA/b_ptcl3D.txt
  echo "stage A dump: $(wc -l < $DA/a_ptcl3D.txt) rows, node $(wc -l < $D/n_ptcl3D.txt)"
  D=$DA
fi
O=$L/runs/$KEY/f2r_$TAG; mkdir -p $O; cd $O
/usr/bin/time -f "wall %e s, maxrss %M kB" $PY $L/f2r_calib.py --node $NODE --dump $D --out $O --cache $L/cache/f2r_${KEY}_${TAG}.npy \
  --truth $TRUTH --nproc $SLURM_CPUS_PER_TASK --cmp $L/dumps/$KEY/b_ptcl3D.txt $EXTRA > $O/f2r.log 2> $O/f2r_time.txt
echo "exit $? $(tail -1 $O/f2r_time.txt)"; tail -40 $O/f2r.log
