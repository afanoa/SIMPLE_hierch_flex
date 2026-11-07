#!/bin/bash
#SBATCH -J klem_lsamp_rec
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=120G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape samples -> full-box maps: plain reconstruct3D (no focus) of lsamp.py's samples at the node's refined poses, lp-10 unsharpened view
# maps (make_view.py, box 160; full maps in VIEW, halves in VIEW/halves) and, for real RNase P (no GT), the closed-/open-arm readout + old
# f692r closed states 15+16 traced into every sample (node_readout2.sbatch) + the radial extent of the closed arm in the end samples.
# env PROJ (project with the node's poses, read-only: copied), LAB (lsamp OUT/samples.txt, master rows), OUT (work dir), VIEW (map folder),
# NAME, MSKDIAM (232), GT (optional: synthetic, maps only)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
O=${OUT:?}; V=${VIEW:?}; NAME=${NAME:-lsamp}; mkdir -p $O $V && cd $O; log(){ echo "[$(date '+%F %T') lsamp_rec $NAME] $*"; }
NS=$(awk '$1 > 0' ${LAB:?} | sort -n | tail -1); log "$NS samples from $LAB, poses $PROJ"
if [ -z "$(ls $O/*_reconstruct3D/recvol_state$(printf %02d $NS)_even_unfil.mrc 2>/dev/null)" ]; then cp -f ${PROJ:?} $O/r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$LAB mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=3 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=${MSKDIAM:-232} nstates=$NS nparts=3 nthr=8 > rec.log 2>&1 || { log "reconstruct FAILED: $(grep -a -m2 ERROR rec.log)"; exit 1; }; fi
RR=$(ls -td $O/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*; : > $O/spec.txt
for s in $(seq 1 $NS); do b=$RR/recvol_state$(printf %02d $s); [ -s ${b}_even_unfil.mrc ] && echo "sample_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec.txt; done
if [ -n "${GT:-}" ]; then $PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1
else RD=$RR LAB=$LAB OUT=$O/readout VIEW=$V NAME=$NAME bash $T/node_readout2.sbatch 2>&1 | grep -v Warning
  CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); TOP=$(sort -k3,3gr $O/readout/states.txt | head -3 | awk '{printf "sample_%02d:%s/state_%02d.mrc ", $1, V, $1}' V=$V)
  $PY $T/radial_extent.py $T/maps160/toy_consensus.mrc $CL --maps $TOP 2>&1 | grep -v "Warning\|lp = lambda" | sed -n '1,8p'; fi
mkdir -p $V/halves; ls $V/*_even.mrc $V/*_odd.mrc > /dev/null 2>&1 && mv -f $V/*_even.mrc $V/*_odd.mrc $V/halves/
log "DONE: full maps in $V ($(ls $V/*.mrc 2>/dev/null | wc -l)), halves in $V/halves"
