#!/bin/bash
#SBATCH -J klem_coord
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# reconstruct the sample cells of the sampled nodes of a coord_tree run (real data) and read them out: per cell the closed-arm /
# open-arm occupancy, particles and its mean position on each axis (coord.txt, coord_axisN.txt). Maps: NODE/maps/state_NN.mrc (full),
# NODE/maps/halves/. env RTAG, FAM, NODESDIR (nodes2), MSKDIAM (232), SNAP
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
O=$K_/${RTAG:?}/coordtree_fam${FAM:?}; N=$O/${NODESDIR:-nodes2}; MD=${MSKDIAM:-232}; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
log(){ echo "[$(date '+%F %T') recon-nodes $RTAG fam$FAM] $*"; }
for nd in $(cat $N/sampled_nodes.txt); do D=$N/node$(printf %03d $nd); [ -s $D/arm_states.txt ] && continue; mkdir -p $D/rec && cd $D/rec
  NB=$(awk '$1 > 0' $D/even_bins.txt | sort -n | uniq | wc -l); cp -f $O/rec/r.simple b.simple; log "node $nd: $NB cells"
  simple_exec prg=selection projfile=b.simple oritype=ptcl3D infile=$D/even_bins.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=b.simple qsys_name=local nparts=4 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=b.simple pgrp=c1 mskdiam=$MD nstates=$NB nparts=4 nthr=8 > rec.log 2>&1 || { log "node $nd reconstruct FAILED"; continue; }
  BR=$(ls -td $D/rec/*_reconstruct3D | head -1); mkdir -p $D/maps/halves
  for s in $(seq 1 $NB); do b=recvol_state$(printf %02d $s); t=state_$(printf %02d $s); cp -f $BR/$b.mrc $D/maps/$t.mrc; cp -f $BR/${b}_even.mrc $D/maps/halves/${t}_even.mrc; cp -f $BR/${b}_odd.mrc $D/maps/halves/${t}_odd.mrc; done
  $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $D/arm_states.txt $D/maps 2>&1 | grep -v Warning | tail -1
  $PY - <<PYEOF
import numpy as np, glob, re
lab = np.loadtxt('$D/even_bins.txt', dtype=int); ax = [np.loadtxt('$D/coord.txt')] + [np.loadtxt(f) for f in sorted(glob.glob('$D/coord_axis*.txt'))]
arm = {int(re.findall(r'state_(\d+)', l.split()[0])[0]): (float(l.split()[1]), float(l.split()[2])) for l in open('$D/arm_states.txt') if l.strip() and not l.startswith('#')}
sd = [np.nanstd(x[lab > 0]) for x in ax]
print('node $nd cells (cell: closed-arm open-arm particles | mean position on axes 1..%d in sd units), by closed-arm:' % len(ax))
for c in sorted(arm, key=lambda c: -arm[c][0]):
    m = lab == c; print('   %02d: %.2f %.2f %6d | %s' % (c, arm[c][0], arm[c][1], m.sum(), ' '.join('%+.2f' % (np.nanmean(x[m]) / s) for x, s in zip(ax, sd))))
PYEOF
done; log DONE
