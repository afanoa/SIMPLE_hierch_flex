#!/bin/bash
#SBATCH -J klem_f692top
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# full692 factorial_v3 (port): top-N maps of chosen site poles ranked by their posterior (POLES 'site:pole:level ..'), plain reconstruct3D
# at the level node's input poses, 10 A views, closed / open arm readout; read-only on the port tree. env POLES N
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K=$A/rnasep_klem; P=$K/stree_port/full692/1_hetero_tree; F=$P/factorial_v3
O=$K/factorial/full692_sites_check/top; V=$K/view_rflex/full692_site_tops; mkdir -p $O $V; : > $O/spec.txt; N=${N:-2500}
for e in ${POLES:?}; do IFS=: read s j lev <<< "$e"; D=$O/${s}_p${j}_top$N; mkdir -p $D
  $PY -c "
import numpy as np
w = [l.split() for l in open('$F/table.txt') if l.startswith('$s ')][0]; Q = np.loadtxt(w[3], ndmin=2); NM = len(Q) if Q.shape[0] > 1000 else Q.shape[1]; Q = Q if Q.shape[0] == NM else Q.T
P = Q[:, $j - 1]; sel = np.zeros(len(P), int); sel[np.argsort(-np.where(P >= 0, P, -np.inf))[:$N]] = 1; np.savetxt('$D/rows.txt', sel, fmt='%d'); print('$s pole $j: top $N by posterior, min P %.3f' % P[sel > 0].min())"
  if [ -z "$(ls $D/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cd $D; cp -f $P/node_$lev/n.simple r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$D/rows.txt mkdir=no > sel.log 2>&1 \
    && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=2 nthr=8 mkdir=no > up.log 2>&1 && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=1 nparts=2 nthr=8 > rec.log 2>&1 \
    || { echo "$s pole $j reconstruct FAILED: $(grep -a -m2 ERROR rec.log)"; cd $O; continue; }; cd $O; fi
  RD=$(ls -td $D/*_reconstruct3D | head -1); echo "${s}_p${j}_top$N $RD/recvol_state01_even_unfil.mrc $RD/recvol_state01_odd_unfil.mrc" >> $O/spec.txt; done
$PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm.txt $(awk -v v=$V '{print v "/" $1 ".mrc"}' $O/spec.txt) 2>&1 | grep -v Warning | tail -1; cat $O/arm.txt
