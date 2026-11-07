#!/bin/bash
#SBATCH -J rnasep_subcls_prior_rec
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Does the 2D expansion sub-class prior improve the docked set? The top N1 particles by (own focal score + shrunk
# sub-class mean) reconstructed at the family poses and scored like the top N1 by own score (same vis.py table).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; FD=${FD:-$KF/closed_focal}
D=$FD/box_subcls/rec; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') sp] $*"; }
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$FD/box_subcls/states_subcls.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for h in "" _even _odd; do cp -f $RD/recvol_state01$h.mrc $D/maps/ 2>/dev/null; done
$PY $FD/lever/vis.py $FD/lever/maps/recvol_state01.mrc "top 30000 by own focal score" $D/maps/recvol_state01.mrc "top 30000 by own + sub-class prior"
log DONE
