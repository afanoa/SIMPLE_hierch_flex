#!/bin/bash
#SBATCH -J rnasep_family_pipeline
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The toy's winning recipe on one real-data family: child board (expand dir) -> k-means K on its posterior means ->
# refine3D_states K states (local 20->8 A) -> refined state maps scored against the deposited conformers -> child board
# again at the refined poses (focus, no free comps, pose marg + devw + pnoise) -> identification table.
# env: E (expansion dir), FAM (1), K (10), MAXITS (6), NP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=${E:-$A/rnasep_klem/expand1}; FAM=${FAM:-1}; K=${K:-10}; NP=${NP:-10}
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; META=$A/rnasep_coupled/basis_span_all/meta_all.npz; CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt
STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
D=$E/fam${FAM}_pipeline; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') fp$FAM] $*"; }
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$E/fam$FAM/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); z = L[:, 1:3]; gm = L[:, 5]; ok = gm > 0.5; rng = np.random.default_rng(0); X = z[ok]; K = $K; best = None
for seed in range(10):
    C = X[rng.choice(len(X), K, replace=False)]
    for it in range(100):
        d = ((X[:, None, :] - C[None]) ** 2).sum(2); a = d.argmin(1); Cn = np.array([X[a == k].mean(0) if (a == k).any() else C[k] for k in range(K)])
        if np.allclose(Cn, C): break
        C = Cn
    cost = d[np.arange(len(X)), a].sum()
    if best is None or cost < best[0]: best = (cost, a.copy())
N = len(np.load('$META')['flags']); st = np.zeros(N, int); st[pind[ok] - 1] = best[1] + 1; np.savetxt('$D/states.txt', st, fmt='%d'); print('family $FAM: k-means K=%d on the child board posterior means: sizes %s' % (K, np.bincount(best[1])))
PYEOF
cp -f $SRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
export OMP_NUM_THREADS=8; log "refine3D_states nstates=$K local 20->8 A maxits ${MAXITS:-6}"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-6} mskdiam=220 pgrp=c1 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_refine3D_states | head -1); mkdir -p $D/maps; for w in $(seq -f %02g 1 $K); do cp -f $RD/recvol_state$w.mrc $D/maps/ 2>/dev/null; done
simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/refined_ptcl3D.txt
$PY - <<PYEOF
import numpy as np, glob
from scipy.ndimage import gaussian_filter
SMPD = 1.464; W = '$W'
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); cons = read_mrc(W + '/maps160/toy_consensus.mrc'); lp = lambda v: gaussian_filter(v, 2.5); L18, Lc = lp(c18), lp(cons)
dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
n = 160; zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); msk = 1. / (1. + np.exp((r - 100 / SMPD) / 2.)); F = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))); band = (r >= n * SMPD / 60.) & (r <= n * SMPD / 20.)
ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; T = {i: F(read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) - cons)[band] for i in ids}
def bc(x, y): return float(np.real(np.sum(np.conj(x) * y)) / np.sqrt(np.sum(np.abs(x) ** 2) * np.sum(np.abs(y) ** 2)))
st = []
for l in open('$D/refined_ptcl3D.txt'):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t); st.append(int(float(kv.get('state', 0))))
st = np.array(st); cover = {}
print('REFINED STATES of family $FAM vs the deposited conformers (full box, 60-20 A):')
for w in range(1, $K + 1):
    V = read_mrc('$D/maps/recvol_state%02d.mrc' % w); Lv = lp(V); occ = (Lv[arm].mean() / Lv[core].mean()) / ref; Dw = F(V - cons)[band]; cc = {i: bc(Dw, T[i]) for i in ids}; best = max(cc, key=cc.get); cover.setdefault(best, []).append((w, cc[best]))
    print('  state %2d: %6d particles | best %d (%+.2f) | closed-arm occ %.2f | %s' % (w, int((st == w).sum()), best, cc[best], occ, ' '.join('%d:%+.2f' % (i, cc[i]) for i in ids)))
print('  conformers identified by a refined state (best match, corr >= 0.3): %s' % ' '.join('%d(%d states, max %+.2f)' % (i, len(v), max(c for _, c in v)) for i, v in sorted(cover.items()) if max(c for _, c in v) >= 0.3))
PYEOF
# feedback: the child board again at the refined poses
$PY - <<PYEOF
import numpy as np
M = dict(np.load('$META')); keys = ('e1', 'e2', 'e3', 'x', 'y'); new = {k: M[k].copy() for k in keys}; n = 0
for i, l in enumerate(open('$D/refined_ptcl3D.txt')):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t)
    if float(kv.get('state', 0)) <= 0: continue
    for k in keys: new[k][i] = float(kv[k])
    n += 1
for k in keys: M[k] = new[k]
np.savez('$D/meta_fb.npz', **M); print('refined poses for %d particles -> meta_fb.npz' % n)
PYEOF
export OMP_NUM_THREADS=1; mkdir -p $D/klem_fb; log "child board at the refined poses"
$PY $W/klem.py $STK none none $D/klem_fb --meta $D/meta_fb.npz --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 0 --nproc 16 --init npz --nsub 0 --subset $E/family${FAM}_subset.txt --ctf auto \
    --lps "14 14 12 12 10 10 8 8 8" --sigmas "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" --focus 1 --pose_marg 3 --devw 0.05 --pnoise --tag "rnasep fam$FAM feedback" 2>&1 | grep --line-buffered -v "^   [0-9]* / " | grep "subset\|final:\|states: cross-half\|Error\|Traceback"
echo "== identification of the feedback child board"; $PY $W/deposited_match.py $D/klem_fb/cells $W/maps160/docked $W/maps160/toy_consensus.mrc
log DONE
