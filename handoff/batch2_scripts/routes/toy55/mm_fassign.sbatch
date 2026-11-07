#!/bin/bash
#SBATCH -J mm_fassign
#SBATCH -p csbdevel
#SBATCH --constraint=x6342
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
# factorial re-assignment of a hetero_tree run (the post-hoc factor layer, tested): every particle scored against the 2^r combination maps
# that mm_leaf_factors.py synth=yes made from the tree's own leaf maps (FDIR/comb_NN.mrc, FDIR/comb_states.txt = initial states), in ONE
# refine3D pass with the tree's stage-B keys (prob_neigh geom, lp 10, box 88, nspace 5000, pcg focus in the root node mask). ITS (1) passes.
# env MSKDIAM (232), PGRP (c1), SMPDCROP (160 x 1.464 / 88: box x smpd / 88). Writes OUT/f.simple, OUT/states.txt. usage: sbatch --export=ALL,RUN=..,FDIR=..,OUT=..,BUILD=.. mm_fassign.sbatch
set -euo pipefail
: ${RUN:?} ${FDIR:?} ${OUT:?} ${BUILD:?}
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
export SIMPLE_PATH=$BUILD PATH=$BUILD/bin:$BUILD/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel SIMPLE_EMAIL=nil OMP_NUM_THREADS=8 SBATCH_CONSTRAINT=x6342
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55
test ! -e $OUT || { echo "$OUT exists"; exit 4; }
mkdir -p $OUT; cd $OUT; cp "$0" launcher.sbatch
NS=$(ls $FDIR/comb_*.mrc | wc -l); MASK=$RUN/1_hetero_tree/node_r/mask.mrc; test -s $MASK
simple_exec prg=print_project_field projfile=$RUN/input.simple oritype=ptcl3D mkdir=no 2>/dev/null | $PY $W/dump_ptcl.py in_dump.npz
$PY $W/set_states_binoris.py $RUN/input.simple f.simple $FDIR/comb_states.txt in_dump.npz
V=""; for j in $(seq 1 $NS); do V="$V vol$j=$FDIR/comb_$(printf %02d $j).mrc"; done
printf '%s\n' "run=$RUN" "fdir=$FDIR" "build=$BUILD" "nstates=$NS" "its=${ITS:-1}" "job=${SLURM_JOB_ID:-}" "started=$(date -Is)" > manifest.txt
simple_exec prg=refine3D projfile=f.simple nstates=$NS $V refine=prob_neigh prob_neigh_mode=geom prob_athres=0.01 lp=10 maxits=${ITS:-1} \
  center=no rec_backend=pcg pcg_mskfile=$MASK pcg_focus=yes nspace=5000 box_crop=88 smpd_crop=${SMPDCROP:-$(echo "160*1.464/88" | bc -l | cut -c1-7)} \
  pgrp=${PGRP:-c1} mskdiam=${MSKDIAM:-232} nparts=4 nthr=8 > refine.log 2>&1 || { echo "refine3D FAILED" >> manifest.txt; exit 7; }
RD=$(ls -td $OUT/*_refine3D 2>/dev/null | head -1); P=${RD:+$RD/}f.simple
simple_exec prg=print_project_field projfile=$P oritype=ptcl3D mkdir=no 2>/dev/null | $PY $W/dump_ptcl.py out_dump.npz
$PY -c "import numpy as np; np.savetxt('states.txt', np.nan_to_num(np.load('out_dump.npz')['state']).astype(int), fmt='%d')"
echo "finished=$(date -Is) states=$(wc -l < states.txt)" >> manifest.txt
