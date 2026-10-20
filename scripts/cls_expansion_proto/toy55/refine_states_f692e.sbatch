#!/bin/bash
#SBATCH -J klem_refine4
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# How far does refinement take the states? 4 states of f692e with very different populations (~39k, ~15k, ~5k leaf states +
# the pure closed zoomed sub-state 4.9k; hard labels, the stored-weights build is not in this snapshot), refine3D_states from the consensus poses against these pure references, local 20 -> LPSTOP A,
# MAXITS iterations, 20 parts x 8 threads; per-state FSC before (the consolidated reconstructions) and after.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692e; O=$R/refine4; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') refine4] $*"; }
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; K=4
$PY - <<PYEOF
import numpy as np
st = np.loadtxt('$R/consolidated_12_fullw_region/states.txt', dtype=int); zs = np.loadtxt('$R/consolidated_12_fullw_zoom/states.txt', dtype=int)
sizes = np.bincount(st); order = [s for s in np.argsort(-sizes) if s > 0]
pick = [order[0], min(order, key=lambda s: abs(sizes[s] - 15000)), min(order, key=lambda s: abs(sizes[s] - 5000))]   # ~39k, ~15k, ~5k particles
out = np.zeros(len(st), int)
for r, s in enumerate(pick): out[st == s] = r + 1
out[zs == 1] = $K   # the pure closed zoomed sub-state (occupancy 0.94, 4.9k) overrides
np.savetxt('$O/states.txt', out, fmt='%d'); print('states: ' + ' '.join('%d:%d' % (r + 1, (out == r + 1).sum()) for r in range($K)) + ' (leaf states %s + closed zoom state)' % pick)
PYEOF
cp -f $SRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=86000 nparts=${NP:-20} nthr=8 mkdir=no > up.log 2>&1
log "refine3D_states nstates=$K pose_policy=local lpstart=20 lpstop=${LPSTOP:-5} maxits=${MAXITS:-10} nparts ${NP:-20} x 8"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=20 lpstop=${LPSTOP:-6} maxits=${MAXITS:-8} mskdiam=220 pgrp=c1 nparts=${NP:-20} nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
grep -a "RESOLUTION @ FSC=0.143" refine.log | tail -3 | cut -c1-160
RD=$(ls -td $O/*_refine3D_states | head -1); mkdir -p $O/maps; cp -f $RD/recvol_state*.mrc $O/maps/ 2>/dev/null; ls $O/maps | head -40 | tr "\n" " "; echo
$PY - <<PYEOF
import numpy as np, glob, os, re, sys
sys.path.insert(0, '$W'); from merge_states_lib import leaf_ok
from model_class import read_mrc
def fsc_res(fe, fo, smpd=1.464):
    e, o = read_mrc(fe).astype(np.float32), read_mrc(fo).astype(np.float32); n = e.shape[0]; x = np.arange(n) - n // 2
    rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
    E, O = np.fft.fftn(e * msk), np.fft.fftn(o * msk); f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), n // 2)
    num = np.bincount(sh.ravel(), np.real(np.conj(E) * O).ravel()); d1 = np.bincount(sh.ravel(), np.abs(E).ravel() ** 2); d2 = np.bincount(sh.ravel(), np.abs(O).ravel() ** 2)
    fsc = num / np.maximum(np.sqrt(d1 * d2), 1e-30); below = np.where(fsc[1:n // 2] < 0.143)[0]; k = (below[0] + 1) if len(below) else n // 2 - 1; return n * smpd / k
st = np.loadtxt('$O/states.txt', dtype=int)
print('state  particles  FSC_after(A)')
for s in range(1, $K + 1):
    fe, fo = '$O/maps/recvol_state%02d_even.mrc' % s, '$O/maps/recvol_state%02d_odd.mrc' % s
    print('%2d  %7d  %s' % (s, (st == s).sum(), '%.1f' % fsc_res(fe, fo) if os.path.isfile(fe) and os.path.isfile(fo) else 'no half maps written'))
PYEOF
CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1); echo "== closed-arm occupancy after refinement"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $O/arm_states.txt $O/maps 2>&1 | grep -v Warning | grep -v mask | grep "state04\|state01\|summary"
log DONE
