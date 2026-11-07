#!/bin/bash
#SBATCH -J mrw_site
#SBATCH -p norm
#SBATCH --constraint=x6740P|x6342
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%A_%a.out
# rerun saved site_rank jobs with another build: LIST has lines "<source site dir> <output dir>"; BUILD=<install>
set -uo pipefail
: ${LIST:?} ${BUILD:?}
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$BUILD/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
export SIMPLE_PATH=$BUILD OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
read -r SRC OUT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)
test -e $SRC/hetero_step.sh || { echo "no hetero_step.sh in $SRC"; exit 3; }
test ! -e $OUT || { echo "$OUT exists"; exit 4; }
mkdir -p $OUT && cd $OUT && cp $SRC/s.simple . || exit 5
cmd=$(grep "prg=site_rank" $SRC/hetero_step.sh | sed "s#^[^ ]*simple_exec#$BUILD/bin/simple_exec#; s#nthr=[0-9.]*#nthr=$SLURM_CPUS_PER_TASK#; s# > hetero_step.log 2>&1##")
printf '%s\n' "src=$SRC" "build=$BUILD" "job=$SLURM_JOB_ID" "cmd=$cmd" "started=$(date -Is)" > rerun_manifest.txt
/usr/bin/time -f "wall %e s, maxrss %M kB" bash -c "$cmd" > site_rank.log 2> time.txt; echo $? > site_rank.exit
echo "finished=$(date -Is) exit=$(cat site_rank.exit) $(tail -1 time.txt)" >> rerun_manifest.txt
