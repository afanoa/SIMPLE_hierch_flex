#!/bin/bash
#SBATCH -J rnasep_level2_pc1
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Truth-free level 2 with the node's OWN main deviation axis: the within-state boards of the two docked-like states
# found without any deposited mask (region_1 state 2, region_p1 state 3) ranked along pc1 of their own cell deviations,
# both ends (10k each) reconstructed at the family poses and scored for the contact against the deposited-mask result.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; KN=$A/rnasep_klem/knobfree; META=$P/meta_fb.npz; VIS=$KF/closed_focal/lever/vis.py
NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))"); PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; log(){ echo "[$(date '+%F %T') l2] $*"; }
for L in $KN/region_1/level2_s2 $KN/region_p1/level2_s3; do
  tag=$(basename $(dirname $L))_$(basename $L); O=$L/pc1; mkdir -p $O; log "$tag: ranking along the own main axis"
  $PY $W/rank_by_target.py $L/board --target pc1 --n 10000 --nrows $NROWS --out $O 2>&1 | grep -v Warning
  mkdir -p $O/rec && cd $O/rec && cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/states_top.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; cd $KN; continue; }
  simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=2 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED"; cd $KN; continue; }
  RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $O/maps; for s in 01 02; do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $O/maps/ 2>/dev/null; done; done; cd $KN
  $PY $VIS $KF/closed_focal_docked30k/rank/maps/recvol_state01.mrc "deposited-mask within-10k (reference)" $L/maps/recvol_state01.mrc "$tag parent-direction top 10k" $O/maps/recvol_state01.mrc "$tag own-axis plus end 10k" $O/maps/recvol_state02.mrc "$tag own-axis minus end 10k" 2>&1 | grep -v Warning
done
log DONE
