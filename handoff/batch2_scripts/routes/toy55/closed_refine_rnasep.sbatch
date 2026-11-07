#!/bin/bash
#SBATCH -J rnasep_closed_refine
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Per-state refinement of the CLOSED state only: the particles whose argmax cell matches conformer 18 (40k), poses from
# the family's 10-state refinement, refined locally against the closed map itself (our conformer-18 map as the start),
# wide mask (240 A), 12 -> 6 A; then the arm occupancy / FSC of the result vs the plain reconstruction (0.76).
# env: LPSTOP (6), MAXITS (6), NP (10), CONF (18)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; MR=$P/klem_fb/matched_rec; CONF=${CONF:-18}
D=$MR/closed_refine; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') cr] $*"; }
$PY - <<PYEOF
import numpy as np
ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; st = np.loadtxt('$MR/states.txt', dtype=int); k = ids.index($CONF) + 1; sel = (st == k).astype(int); np.savetxt('$D/states.txt', sel, fmt='%d'); print('closed particles (argmax cell matched to %d): %d' % ($CONF, sel.sum()))
PYEOF
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
K=$(printf %02d $(python3 -c "print([1,10,18,27,28,37,46,48,49,51].index($CONF)+1)")); cp -f $MR/maps/recvol_state$K.mrc $D/start_closed.mrc
log "refine3D single state, local poses against the closed map, mskdiam 240, 12 -> ${LPSTOP:-6} A, maxits ${MAXITS:-6}"
simple_exec prg=refine3D projfile=p.simple vol1=$D/start_closed.mrc pose_policy=local lpstart=12 lpstop=${LPSTOP:-6} maxits=${MAXITS:-6} mskdiam=240 pgrp=c1 nparts=${NP:-10} nthr=8 > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; tail -5 refine.log; exit 1; }
RD=$(ls -td $D/*_refine3D | head -1); ls $RD | grep -i 'recvol\|rec_final\|pproc' | head -8; mkdir -p $D/maps; for h in "" _even _odd; do cp -f $RD/recvol_state01$h.mrc $D/maps/ 2>/dev/null; done; cp -f $RD/rec_final_state01*.mrc $D/maps/ 2>/dev/null
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464; S = Shells(160, SMPD); sh = S.sh; msk = S.soft_mask(116 / SMPD); g = lambda v: gaussian_filter(v, 2.5)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
def fsc_res(x, y):
    a_, b_ = fftc(x * msk), fftc(y * msk); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return 160 * SMPD / s
for name in sorted(glob.glob('$D/maps/*.mrc')):
    if '_even' in name or '_odd' in name: continue
    v = read_mrc(name); occ = (g(v)[arm].mean() / g(v)[core].mean()) / ref
    try: res = fsc_res(read_mrc(name.replace('.mrc', '_even.mrc')), read_mrc(name.replace('.mrc', '_odd.mrc')))
    except Exception: res = float('nan')
    print('CLOSED REFINED %s: closed-arm occupancy %.2f | FSC 0.143 %.1f A' % (name.split('/')[-1], occ, res))
v = read_mrc('$D/start_closed.mrc'); print('start map (plain reconstruction): closed-arm occupancy %.2f' % ((g(v)[arm].mean() / g(v)[core].mean()) / ref))
PYEOF
RDP=$RD; simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/refined_ptcl3D.txt; log "refined table rows: $(grep -c . $D/refined_ptcl3D.txt)"
log DONE
