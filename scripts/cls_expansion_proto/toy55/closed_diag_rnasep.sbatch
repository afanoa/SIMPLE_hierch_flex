#!/bin/bash
#SBATCH -J rnasep_closed_diag
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Why is the closed arm short? (a) arm extent vs the masks used (100 A soft mask in scoring/sharpening, mskdiam 220 in
# reconstruct3D); (b) per-cell arm occupancy on family 1's feedback board; (c) nested reconstructions of the closed
# particles ranked by their responsibility for the 18-matched cells (top 10k, next 10k, next 20k) at the refined poses:
# does the arm grow with confidence (mixture) or not (blur / disorder)?
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KD=$P/klem_fb; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
D=$KD/closed_diag; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') cd] $*"; }
echo "== (a) arm extent"; $PY - <<PYEOF
import numpy as np, glob
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); n = 160; g = lambda v: gaussian_filter(v, 2.5)
L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18
zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * SMPD
print('arm voxels %d; radius from the box centre: median %.0f A, 90th pct %.0f A, max %.0f A | dense-18 voxels max radius %.0f A | scoring/sharpening soft mask radius 100 A, reconstruct3D mskdiam 220 A (radius 110)' % (arm.sum(), np.median(r[arm]), np.percentile(r[arm], 90), r[arm].max(), r[dens18].max()))
print('fraction of arm voxels beyond 95 A: %.2f, beyond 100 A: %.2f, beyond 110 A: %.2f' % (np.mean(r[arm] > 95), np.mean(r[arm] > 100), np.mean(r[arm] > 110)))
PYEOF
echo "== (b) per-cell closed-arm occupancy on the feedback board"; $PY $W/closed_arm_cells.py $KD --top 8
echo "== (c) nested reconstructions of the closed particles by confidence"; $PY - <<PYEOF
import numpy as np
ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; M = np.loadtxt('$KD/matched_rec/cell_matches.txt', comments='#'); best = {int(r[0]): int(r[1]) for r in M}
R = np.load('$KD/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); ngrid = len(R['i1']); cells18 = [k for k in range(ngrid) if best.get(k + 1) == 18]; cellsC = [k for k in range(ngrid) if best.get(k + 1) in (18, 46, 10)]
s18 = r[:, cells18].sum(1); sC = r[:, cellsC].sum(1); order = np.argsort(-s18); N = len(np.load('$META')['flags']); st = np.zeros(N, int)
for j, (lo, hi) in enumerate(((0, 10000), (10000, 20000), (20000, 40000))): st[pinds[order[lo:hi]] - 1] = j + 1
np.savetxt('$D/states.txt', st, fmt='%d'); print('cells matched to 18: %d, to the closed family (18/46/10): %d; closed-18 responsibility: top 10k >= %.2f, 10-20k >= %.2f, 20-40k >= %.2f (mean %.3f)' % (len(cells18), len(cellsC), s18[order[9999]], s18[order[19999]], s18[order[39999]], s18.mean()))
PYEOF
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=240 nstates=3 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for w in 01 02 03; do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps/ 2>/dev/null; done; done
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464; S = Shells(160, SMPD); sh = S.sh; msk = S.soft_mask(110 / SMPD); g = lambda v: gaussian_filter(v, 2.5)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
def fsc_res(x, y):
    a_, b_ = fftc(x * msk), fftc(y * msk); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return 160 * SMPD / s
V = [read_mrc('$D/maps/recvol_state%02d.mrc' % w) for w in (1, 2, 3)]; nn = [10000, 10000, 20000]; cum = None; ncum = 0
for j, (v, k) in enumerate(zip(V, nn)):
    occ = (g(v)[arm].mean() / g(v)[core].mean()) / ref; cum = v * k if cum is None else cum + v * k; ncum += k; oc = (g(cum / ncum)[arm].mean() / g(cum / ncum)[core].mean()) / ref
    try: res = fsc_res(read_mrc('$D/maps/recvol_state%02d_even.mrc' % (j + 1)), read_mrc('$D/maps/recvol_state%02d_odd.mrc' % (j + 1)))
    except Exception: res = float('nan')
    print('  shell %d (%s): closed-arm occupancy %.2f (FSC %.1f A) | cumulative top %dk: %.2f' % (j + 1, ('top 10k', 'next 10k', 'next 20k')[j], occ, res, ncum // 1000, oc))
PYEOF
log DONE
