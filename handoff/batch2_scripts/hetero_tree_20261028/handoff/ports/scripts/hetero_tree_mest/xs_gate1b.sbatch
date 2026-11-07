#!/bin/bash
#SBATCH -J xs_gate1b
#SBATCH -p csbdevel
#SBATCH --constraint=x6342
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest/logs/%x_%A_%a.out
# site_split gate 1: the worker alone on one saved gate node per array task (LIST lines: "KEY NODE MSKDIAM PGRP"), every site of
# the node's sites/ (and islands when present); summary.txt per node under $M/xs_gate1/KEY for the comparison with the Python calls
set -uo pipefail
: ${LIST:?} ${BUILD:?}
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$BUILD/bin:$BUILD/scripts:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
export SIMPLE_PATH=$BUILD SIMPLE_SPLIT_CONTRAST=maps OMP_NUM_THREADS=16
read -r KEY NODE MSK PG < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
test -s $NODE/stageB/b.simple || { echo "no stage-B project at $NODE"; exit 3; }
W=$M/${OUT:-xs_gate1}/$KEY; test ! -e $W || { echo "$W exists"; exit 4; }; mkdir -p $W; cd $W
cp $NODE/stageB/b.simple x.simple
echo "node $NODE: $(ls $NODE/sites/site_*.mrc 2>/dev/null | wc -l) sites, $(ls $NODE/sites/island_*.mrc 2>/dev/null | wc -l) islands; started $(date -Is)"
/usr/bin/time -f "wall %e s, maxrss %M kB" simple_exec prg=site_split projfile=x.simple dir=$NODE mskdiam=$MSK pgrp=$PG nthr=16 mkdir=no > site_split.log 2> time.txt
echo "exit $? $(tail -1 time.txt)"; grep -h "SITE_SPLIT\|cross-half EM\|renumbered" site_split.log | tail -40; echo "--- summary"; cat summary.txt
