#!/bin/bash
#SBATCH -J klem_states_speed
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# timing of the chain's family refine3D_states step on an existing family input (env SRCDIR = <run>/fam<j>/refine with
# p.simple carrying the state labels, TAG, EXTRA = extra refine3D_states arguments, e.g. "filt_mode=none"); same
# settings as f692_family.sh (nstates from the project, local poses 20 -> 8 A, maxits 6, 10 parts x 8 threads)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; SRCDIR=${SRCDIR:?}; TAG=${TAG:?}; EXTRA=${EXTRA:-}
B=$(grep -h -o "B=\$HOME/SIMPLE_kernel/snap/[a-z0-9]*" $W/f692_family_job3.sbatch $W/f692_family.sh $W/f692_region_one3.sbatch 2>/dev/null | head -1 | sed "s/B=//; s#\$HOME#$HOME#"); B=${B:-$HOME/SIMPLE_kernel/snap/5e1e6aeee}
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
R=$A/rnasep_klem/states_speed/$TAG; mkdir -p $R; cd $R; cp -f $SRCDIR/p.simple p.simple; log(){ echo "[$(date '+%F %T') speed $TAG] $*"; }
K=$(grep -o "nstates=[0-9]*" $SRCDIR/refine.log 2>/dev/null | head -1 | cut -d= -f2); K=${K:-$(grep -o "NSTATES FROM PROJECT: *[0-9]*" $SRCDIR/refine.log | grep -o "[0-9]*$")}
log "build $B, input $SRCDIR/p.simple, nstates $K, extra: '$EXTRA'"; t0=$(date +%s)
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=20 lpstop=8 maxits=6 mskdiam=220 pgrp=c1 nparts=10 nthr=8 $EXTRA > refine.log 2>&1; rc=$?
log "refine3D_states rc $rc in $(( $(date +%s) - t0 )) s"
grep -a "ITERATION\|NORMAL STOP\|ERROR\|FINAL" refine.log | head -12 | cut -c1-140
d=$(ls -td $R/*_refine3D_states | head -1); for f in $d/VOLASSEMBLE_BENCH_ITER000.txt $d/VOLASSEMBLE_BENCH_ITER003.txt; do [ -f $f ] && { echo "-- $(basename $f)"; grep "nonuniform_filter\|restore_eos\|total time" $f; }; done
ls --time-style=+%H:%M:%S -l $d | awk '{print $6, $7}' | grep "VOLASSEMBLE_BENCH_ITER00[16]\|REC3D_BENCH_ITER00[0789]_PART01\|recvol_state01.mrc$\|RESOLUTION_STATE01$" | sort
log DONE
