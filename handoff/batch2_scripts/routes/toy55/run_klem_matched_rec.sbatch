#!/bin/bash
#SBATCH -J rnasep_matched_rec
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The ten deposited conformers recovered from the real data: on a klem board that identifies them (default: family 1's
# feedback board), every cell is labelled by the deposited conformer it matches best; particles inherit their argmax
# cell's conformer; one full-box reconstruct3D per conformer at the REFINED poses (the family's refine3D_states project);
# every map scored against its conformer (closed-arm occupancy, deviation corr at 60-20 A, FSC resolution).
# env: KDIR (board), PROJSRC (project with the poses), NP (10), MSKDIAM (220), RTAG
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline
KDIR=${KDIR:-$P/klem_fb}; PROJSRC=${PROJSRC:-$(ls -td $P/*_refine3D_states | head -1)/p.simple}; META=$A/rnasep_coupled/basis_span_all/meta_all.npz; NP=${NP:-10}
D=$KDIR/matched_rec${RTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') mrec] $*"; }
$PY $W/deposited_match.py $KDIR/cells $W/maps160/docked $W/maps160/toy_consensus.mrc --save $D/cell_matches.txt | tail -3
$PY - <<PYEOF
import numpy as np
ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; M = np.loadtxt('$D/cell_matches.txt', comments='#'); best = {int(r[0]): int(r[1]) for r in M}
L = np.loadtxt('$KDIR/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); gm = L[:, 5]; ok = (gm > 0.5) & (am <= len(best))
N = len(np.load('$META')['flags']); st = np.zeros(N, int)
for p, c in zip(pind[ok], am[ok]): st[p - 1] = ids.index(best[c]) + 1
np.savetxt('$D/states.txt', st, fmt='%d'); cnt = np.bincount(st, minlength=11)
print('states = deposited conformer of the argmax cell: ' + ' '.join('%d:%d' % (ids[k - 1], cnt[k]) for k in range(1, 11)) + ' | cells per conformer: ' + ' '.join('%d:%d' % (i, sum(1 for v in best.values() if v == i)) for i in ids))
PYEOF
cp -f $PROJSRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; tail -2 sel.log; exit 1; }
log "reconstruct3D 10 states at the refined poses ($PROJSRC)"
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=${MSKDIAM:-220} nstates=10 nparts=$NP nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for w in $(seq -f %02g 1 10); do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps/ 2>/dev/null; done; done
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
SMPD = 1.464; ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; S = Shells(160, SMPD); msk = S.soft_mask(100 / SMPD); sh = S.sh
cons = read_mrc('$W/maps160/toy_consensus.mrc'); lp = lambda v: gaussian_filter(v, 2.5)
c18 = read_mrc(glob.glob('$W/maps160/docked/conf18_*.mrc')[0]); L18, Lc = lp(c18), lp(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
T = {i: fftc((read_mrc(glob.glob('$W/maps160/docked/conf%02d_*.mrc' % i)[0]) - cons) * msk) for i in ids}
def bc(x, y, lo, hi):
    m = S.band(lo, hi); return float(np.real(np.sum(np.conj(x[m]) * y[m])) / max(np.sqrt(np.sum(np.abs(x[m]) ** 2) * np.sum(np.abs(y[m]) ** 2)), 1e-30))
def fsc_res(x, y):
    a_, b_ = fftc(x * msk), fftc(y * msk); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return 160 * SMPD / s
st = np.loadtxt('$D/states.txt', dtype=int); print('RECOVERED CONFORMERS (one full-box map per deposited conformer, refined poses, 60-20 A deviation corr vs every deposited conformer):')
for k, i in enumerate(ids):
    V = read_mrc('$D/maps/recvol_state%02d.mrc' % (k + 1)); Lv = lp(V); occ = (Lv[arm].mean() / Lv[core].mean()) / ref; Dw = fftc((V - cons) * msk); cc = {j: bc(Dw, T[j], 60, 20) for j in ids}; own = cc[i]; bestj = max(cc, key=cc.get)
    try: res = fsc_res(read_mrc('$D/maps/recvol_state%02d_even.mrc' % (k + 1)), read_mrc('$D/maps/recvol_state%02d_odd.mrc' % (k + 1)))
    except Exception: res = float('nan')
    print('  conformer %2d: %6d particles | corr with itself %+.2f (best %d %+.2f) | closed-arm occ %.2f | FSC 0.143 %.1f A' % (i, int((st == k + 1).sum()), own, bestj, cc[bestj], occ, res))
PYEOF
log DONE
