#!/bin/bash
#SBATCH -J land_rec
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=80G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_rec.%j.out
# plain SIMPLE reconstructions of coordinate windows (one state per window) at a node's poses; 10 A unsharpened views (halves/),
# arm readout, radial extent, composition in deposited conformers. env: OUT VIEW LAB (state per master row) PROJ NAME
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 SIMPLE_EMAIL=nil
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
O=${OUT:?}; V=${VIEW:?}; NAME=${NAME:-land}; mkdir -p $O $V/halves && cd $O; log(){ echo "[$(date '+%F %T') land_rec $NAME] $*"; }
NS=$(sort -n ${LAB:?} | tail -1); log "$NS windows from $LAB, poses $PROJ"
if [ -z "$(ls $O/*_reconstruct3D/recvol_state$(printf %02d $NS)_odd_unfil.mrc 2>/dev/null)" ]; then cp -f ${PROJ:?} $O/r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$LAB mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=3 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=$NS nparts=3 nthr=8 > rec.log 2>&1 || { log "reconstruct FAILED: $(grep -a -h -m2 ERROR sel.log up.log rec.log)"; exit 1; }; fi
RR=$(ls -td $O/*_reconstruct3D | head -1); : > $O/spec.txt; CNT=""
for s in $(seq 1 $NS); do b=$RR/recvol_state$(printf %02d $s); [ -s ${b}_even_unfil.mrc ] || continue
  echo "state_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec.txt; cp -f ${b}_even_unfil.mrc $V/halves/state_$(printf %02d $s)_even.mrc; cp -f ${b}_odd_unfil.mrc $V/halves/state_$(printf %02d $s)_odd.mrc
  CNT="$CNT,$(awk -v s=$s '$1 == s' $LAB | wc -l)"; done
$PY $T/make_view.py $V $O/spec.txt --lp 10 --box 160 > /dev/null 2>&1; log "views: $(ls $V/state_*.mrc | wc -l) maps lp 10 unsharpened in $V, halves in $V/halves"
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm.txt $V 2>&1 | grep --line-buffered -v "Warning\|lp = lambda" | sed 's#^#  #'
$PY $T/radial_extent.py $T/maps160/toy_consensus.mrc $CL --maps $(for f in $V/state_*.mrc; do echo "$(basename $f .mrc):$f"; done) 2>&1 | grep --line-buffered -v "Warning\|lp = lambda" | sed -n '1,40p'
$PY $T/land_confmix.py $T/maps160/toy_consensus.mrc $K_/landscape/confs/arm_confs.txt $O/confmix.txt --refs $T/maps160/docked_all/conf*_m160.mrc $T/maps160/holo_docked/conf*_m160.mrc \
   --maps $V/state_*.mrc --counts ${CNT#,} 2>&1 | grep --line-buffered -v Warning
log "DONE: $O (arm.txt, confmix.txt), maps $V"
