#!/bin/bash
#SBATCH -J klem_purify
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# purification of a consolidated landscape: the merged states as starting labels, refine3D_states re-assigns every particle
# (prob_neigh, local poses, 10->8 A, maxits 5) on a particle subset (env SUBSET rows file, default all), then the arm readout
# of the output state maps. env: RTAG, CUT (0.7), CTAG, SUBSET, NP (20), MAXITS (5)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=32
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}; C=$R/consolidated_${CUT:-0.7}${CTAG:-}; P=$C/purify${SUBSET:+_sub}; mkdir -p $P; cd $P; log(){ echo "[$(date '+%F %T') purify] $*"; }
SRC=${SRC:-$A/rnasep_coupled/full692k_2d/rnasep.simple}
$PY - <<PYEOF
import numpy as np
st = np.loadtxt('$C/states.txt', dtype=int); sub = '$SUBSET'
if sub: m = np.zeros(len(st), bool); m[np.loadtxt(sub, dtype=int).ravel() - 1] = True; st = np.where(m, st, 0)
sizes = np.bincount(st)[1:]; keep = np.where(sizes >= 500)[0] + 1; remap = np.zeros(st.max() + 1, int); remap[keep] = np.arange(1, len(keep) + 1); st = remap[st]
np.savetxt('$P/states_in.txt', st, fmt='%d'); print('purify input: %d states (>= 500 particles), %d particles: %s' % (len(keep), (st > 0).sum(), ' '.join(str(int(x)) for x in np.bincount(st)[1:])))
PYEOF
K=$(sort -n $P/states_in.txt | tail -1); cp -f $SRC p.simple
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$P/states_in.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
log "refine3D_states nstates=$K local 10->8 A maxits ${MAXITS:-5} nparts ${NP:-20}"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=10 lpstop=8 maxits=${MAXITS:-5} mskdiam=232 pgrp=c1 nparts=${NP:-20} nthr=32 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-160)"; exit 1; }
RD=$(ls -td $P/*_refine3D_states | head -1); mkdir -p $P/maps; cp -f $RD/recvol_state??.mrc $P/maps/ 2>/dev/null; simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); print a[2]}}' > $P/states_out.txt
$PY -c "import numpy as np; a=np.loadtxt('$P/states_in.txt',dtype=int); b=np.loadtxt('$P/states_out.txt',dtype=int); n=min(len(a),len(b)); a,b=a[:n],b[:n]; m=(a>0)&(b>0); print('purify: %d of %d particles changed state; output sizes %s' % ((a[m]!=b[m]).sum(), m.sum(), ' '.join(str(int(x)) for x in np.bincount(b)[1:])))"
echo "== arm readout of the purified states"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $(ls $W/maps160/docked/conf18_*_m160.mrc | head -1) $P/arm_states.txt $P/maps 2>&1 | grep -v Warning
log "DONE $P"
