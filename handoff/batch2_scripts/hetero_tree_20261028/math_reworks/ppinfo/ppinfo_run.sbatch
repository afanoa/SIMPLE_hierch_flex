#!/bin/bash
#SBATCH -J ppinfo
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=28G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo/logs/%x_%A_%a.out
# one ppinfo_site.py configuration per array task. LIST lines: "OUT NODEKEY NODE SITE TRUTH RARE [python args...]"
# (SITE = a site_rank dir with hetero_step.sh, or a replay dir with rerun_manifest.txt)
#   OUT = output dir name under $P/runs, NODEKEY = the dump / plane-cache key of the node, NODE = the node dir (n.simple, stageB/,
#   sites/, fflex/), SITE = the saved site_rank dir (its hetero_step.sh gives pcg_mskfile, site_dirs and mskdiam; site_contrast.txt the
#   contrast weights), TRUTH = fixture truth.txt or '-', RARE = rare conformer or 0. The node's projects are dumped once (print_project_field
#   on the build's simple_exec) and its planes cached once at box CACHEBOX (88) under $P/cache.
set -uo pipefail
: ${LIST:?}
P=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
B=${BUILD:-/mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim2/build_vd/install}
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
read -r OUT KEY NODE SITE TRUTH RARE EXTRA < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST | sed 's/  */ /g' | awk '{o=$1; k=$2; n=$3; s=$4; t=$5; r=$6; $1=$2=$3=$4=$5=$6=""; print o, k, n, s, t, r, $0}')
test -d $NODE || { echo "no node $NODE"; exit 3; }
D=$P/dumps/$KEY; mkdir -p $D $P/cache $P/runs
if [ ! -s $D/n_ptcl3D.txt ]; then
  simple_exec prg=print_project_field projfile=$NODE/n.simple oritype=stk mkdir=no > $D/n_stk.txt 2> $D/n_stk.err
  simple_exec prg=print_project_field projfile=$NODE/n.simple oritype=ptcl3D mkdir=no > $D/n_ptcl3D.txt 2> $D/n_ptcl3D.err
  [ -s $NODE/stageB/b.simple ] && simple_exec prg=print_project_field projfile=$NODE/stageB/b.simple oritype=ptcl3D mkdir=no > $D/b_ptcl3D.txt 2> $D/b_ptcl3D.err
  echo "dumped $KEY: $(wc -l < $D/n_ptcl3D.txt) ptcl3D rows, $(wc -l < $D/n_stk.txt) stacks"
fi
CMD=$(cat $SITE/hetero_step.sh $SITE/rerun_manifest.txt 2>/dev/null | grep "prg=site_rank" | head -1)
MASK=$(echo "$CMD" | grep -o "pcg_mskfile=[^ ]*" | cut -d= -f2)
DIRS=$(echo "$CMD" | grep -o "site_dirs=[^ ]*" | cut -d= -f2)
MSKD=$(echo "$CMD" | grep -o "mskdiam=[^ ]*" | cut -d= -f2)
test -n "$MASK" -a -n "$MSKD" || { echo "no site_rank command with pcg_mskfile and mskdiam in $SITE"; exit 3; }
TR=""; [ "$TRUTH" != "-" ] && TR="--truth $TRUTH --rare $RARE"
CB=$(echo "$EXTRA" | grep -o "CB=[0-9]*" | cut -d= -f2); CB=${CB:-${CACHEBOX:-88}}; EXTRA=$(echo "$EXTRA" | sed "s/CB=[0-9]*//")
O=$P/runs/$OUT; mkdir -p $O
echo "config: out=$OUT key=$KEY node=$NODE site=$SITE mask=$MASK dirs=$DIRS mskdiam=$MSKD truth=$TRUTH rare=$RARE extra=$EXTRA job=$SLURM_JOB_ID"
cd $O
/usr/bin/time -f "wall %e s, maxrss %M kB" $PY $P/ppinfo_site.py --node $NODE --dump $D --site $SITE --mask $MASK --dirs ${DIRS:-1,3} --mskdiam $MSKD \
  --out $O --cache $P/cache/${KEY}_B$CB.npy --cachebox $CB --nproc $SLURM_CPUS_PER_TASK $TR $EXTRA > run.log 2> time.txt
echo "exit $? $(tail -1 time.txt)"; tail -8 run.log
