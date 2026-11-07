#!/bin/bash
#SBATCH -J state_recs
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# full-size maps of a final labelling (both halves; SIMPLE reconstruct3D per half on the half project, states from the
# labels), the two halves averaged per state, then the closed-arm readout (arm_readout.py) of every state map.
# env MODE, TAG, FIN (final dir, final_xem), CLOSED (deposited closed map in our frame), CONSMAP (consensus)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:?}; FIN=${FIN:-final_xem}
O=$K/dtree_${MODE}_$TAG/$FIN; M=$O/maps; mkdir -p $M; log(){ echo "[$(date '+%F %T') state_recs $MODE $TAG] $*"; }
for H in A B; do D=$O/rec_$H; mkdir -p $D; cd $D; ns=$($PY -c "import numpy as np; print(int(np.loadtxt('$O/children_$H.txt').max()))")
  [ -s $D/done ] || { cp -f $K/synth_${MODE}_$H/src.simple p.simple && simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/children_$H.txt mkdir=no > sel.log 2>&1 && \
    simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=1 nthr=16 > rec.log 2>&1 && touch done || { log "rec $H FAILED: $(grep -a -m2 'ERROR\|STOP' rec.log sel.log | cut -c1-160)"; exit 1; }; }
  log "half $H: $ns states reconstructed"; done
RA=$(ls -td $O/rec_A/*_reconstruct3D 2>/dev/null | head -1); RB=$(ls -td $O/rec_B/*_reconstruct3D 2>/dev/null | head -1); [ -z "$RA" ] && RA=$O/rec_A; [ -z "$RB" ] && RB=$O/rec_B
$PY - <<PYEOF
import sys, glob, os, numpy as np
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
fa = sorted(glob.glob('$RA/recvol_state??.mrc')); n = 0
for f in fa:
    g = f.replace('$RA', '$RB')
    if not os.path.isfile(g): continue
    write_mrc('$M/state_%s.mrc' % os.path.basename(f)[12:14], 0.5 * (read_mrc(f) + read_mrc(g)), 1.464); n += 1
print('averaged %d state maps (halves A + B) into $M' % n)
PYEOF
$PY $W/arm_readout.py ${CONSMAP:?} ${CLOSED:?} $O/arm_readout.txt $M 2>&1 | grep -v Warning
log DONE
