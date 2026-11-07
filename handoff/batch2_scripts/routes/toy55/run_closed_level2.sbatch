#!/bin/bash
#SBATCH -J rnasep_closed_level2
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Level 2 of the real-data tree: the particles of refined state STATE (closed-enriched family) at their REFINED poses ->
# klem on that subset (no free comps, best configuration) -> identification table + PC1 bins -> K-state refinement seeded
# by k-means of the posterior means -> closed-arm occupancy and conformer matches per refined sub-state.
# env: L1 (level-1 refine dir), STATE (1), K (3), LPS, NP
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; NP=${NP:-10}; STATE=${STATE:-1}; K=${K:-3}
L1=${L1:-$A/rnasep_klem/closed_refine_nb3}; RD=$(ls -td $L1/*_refine3D_states | head -1); META=$A/rnasep_coupled/basis_span_all/meta_all.npz; CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt
STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
D=$A/rnasep_klem/closed_level2_s$STATE; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') l2] $*"; }
simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/l1_refined_ptcl3D.txt
$PY - <<PYEOF
import numpy as np
M = dict(np.load('$META')); keys = ('e1', 'e2', 'e3', 'x', 'y'); new = {k: M[k].copy() for k in keys}; sub = []; n = 0
for i, l in enumerate(open('$D/l1_refined_ptcl3D.txt')):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t); st = int(float(kv.get('state', 0)))
    if st <= 0: continue
    for k in keys: new[k][i] = float(kv[k])
    n += 1
    if st == $STATE: sub.append(i + 1)
for k in keys: M[k] = new[k]
np.savez('$D/meta_l1.npz', **M); np.savetxt('$D/subset_state$STATE.txt', np.array(sub), fmt='%d'); print('refined poses for %d particles; state $STATE subset: %d particles' % (n, len(sub)))
PYEOF
export OMP_NUM_THREADS=1; mkdir -p $D/klem
$PY $W/klem.py $STK none none $D/klem --meta $D/meta_l1.npz --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 0 --nproc 16 --init npz --nsub 0 --subset $D/subset_state$STATE.txt --ctf auto \
    --lps "${LPS:-20 20 14 14 10 10 8 8 8}" --sigmas "${SIGMAS:-1.5 1.5 1.0 1.0 0.8 0.7 0.7 0.7 0.7}" --pose_marg 3 --devw 0.05 --pnoise --tag "closed level2 s$STATE" 2>&1 | grep --line-buffered -v "^   [0-9]* / "
echo "== identification of the level-2 states"; $PY $W/deposited_match.py $D/klem/cells $W/maps160/docked $W/maps160/toy_consensus.mrc
# K-state refinement seeded by k-means on the level-2 posterior means (the particles outside the subset keep state 0)
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$D/klem/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); z = L[:, 1:3]; gm = L[:, 5]; ok = gm > 0.5; rng = np.random.default_rng(0); X = z[ok]; K = $K; best = None
for seed in range(10):
    C = X[rng.choice(len(X), K, replace=False)]
    for it in range(100):
        d = ((X[:, None, :] - C[None]) ** 2).sum(2); a = d.argmin(1); Cn = np.array([X[a == k].mean(0) if (a == k).any() else C[k] for k in range(K)])
        if np.allclose(Cn, C): break
        C = Cn
    cost = d[np.arange(len(X)), a].sum()
    if best is None or cost < best[0]: best = (cost, a.copy())
N = len(np.load('$META')['flags']); st = np.zeros(N, int); st[pind[ok] - 1] = best[1] + 1; np.savetxt('$D/states_l2.txt', st, fmt='%d'); print('level-2 k-means K=%d: sizes %s' % (K, np.bincount(best[1])))
PYEOF
cp -f $RD/p.simple $D/p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states_l2.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
export OMP_NUM_THREADS=8; log "refine3D_states nstates=$K local 20->8 A on the level-2 partition"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-6} mskdiam=220 pgrp=c1 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
RD2=$(ls -td $D/*_refine3D_states | head -1); mkdir -p $D/maps; for w in $(seq -f %02g 1 $K); do cp -f $RD2/recvol_state$w.mrc $D/maps/ 2>/dev/null; done
simple_exec prg=print_project_field projfile=$RD2/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); print int(a[2])}}' | sort -n | uniq -c | awk '{printf "state %d: %d particles\n", $2, $1}'
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
M = [read_mrc('$D/maps/recvol_state%02d.mrc' % w) for w in range(1, $K + 1)]
def bc(x, y): return float(np.real(np.sum(np.conj(x) * y)) / np.sqrt(np.sum(np.abs(x) ** 2) * np.sum(np.abs(y) ** 2)))
for w, V in enumerate(M):
    Lv = lp(V); occ = (Lv[arm].mean() / Lv[core].mean()) / ref; Dw = F(V - cons)[band]; cc = {i: bc(Dw, T[i]) for i in ids}; best = max(cc, key=cc.get)
    print('level-2 refined state %d: closed-arm occupancy %.2f (level-1 state 0.66, consensus 0.50, deposited 18 = 1.00) | best conformer %d (%+.2f) | all: %s' % (w + 1, occ, best, cc[best], ' '.join('%d:%+.2f' % (i, cc[i]) for i in ids)))
PYEOF
log DONE
