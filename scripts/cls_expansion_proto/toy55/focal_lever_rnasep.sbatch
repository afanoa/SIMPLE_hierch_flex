#!/bin/bash
#SBATCH -J rnasep_focal_lever
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Closed state from the family-1-wide focal board (409k particles): (1) every particle ranked by its responsibility-
# weighted focal score (sum_m r_im x (corr closed - corr open)_m), the docked population estimated from the score
# distribution, the top N1 / next N2 / bottom N1 reconstructed at the family poses (unmasked) and scored (arm, tip,
# base, FSC, contact visibility at 8 A); (2) LEVER-ARM pose refinement of the top N1 as ONE state: refine3D
# refine=neigh (neighbourhood of ATHRES deg around the family pose, no global search, no re-classification) against the
# docked set's own reconstruction under the full 232 A mask (the arm at 95-115 A radius INCLUDED, the opposite of the
# core-masked run), 10 -> LPSTOP A, gold-standard halves; unmasked reconstruction at the refined poses, scored the same
# way, pose change reported, the refined docked map sharpened (B -150) for download.
# env: FD (board dir), N1 (30000), N2 (20000), ATHRES (8), LPSTOP (6), MAXITS (6), NP (6), RTAG
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; FD=${FD:-$KF/closed_focal}; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
N1=${N1:-30000}; N2=${N2:-20000}; D=$FD/lever${RTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') lv] $*"; }
SC=$KF/closed_focal_closedfam/rank/score.py
cat > $D/vis.py <<'PYEOF'
import numpy as np, glob, sys
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, write_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
W = '/mnt/beegfs/elmlund/afan/routes/toy55'; SMPD = 1.464; n = 160; S = Shells(n, SMPD); sh = S.sh; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); rA = rr * SMPD; s_inv = rr / (n * SMPD)
lp = lambda v, res: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))) * (rr <= n * SMPD / res)))))
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); c27 = read_mrc(glob.glob(W + '/maps160/docked/conf27_*.mrc')[0])
g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); tip = arm & (rA > 95); ref = L18[arm].mean() / L18[core].mean()
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
def row(name, vv, ev):
    L = lp(vv, 8.); L = L / max(np.median(L[core]), 1e-30); Lv = g(vv); occ = (Lv[arm].mean() / Lv[core].mean()) / ref
    try: res = fsc_res(read_mrc(ev[0]), read_mrc(ev[1])) if ev else float('nan')
    except Exception: res = float('nan')
    print('  %-50s %8.2f %8.2f %9.2f %6.2f %6.1f' % (name, np.percentile(L[tip], 99), np.percentile(L[tip], 90), np.mean(L[tip] > 0.5), occ, res))
print('contact VISIBILITY at 8 A (core median = 1) and arm occupancy:  %-26s %8s %8s %9s %6s %6s' % ('map', 'tip p99', 'tip p90', 'tip>0.5', 'arm', 'FSC'))
row('deposited 18 (closed)', c18, None); row('deposited 27 (open)', c27, None)
for p_, nm in zip(sys.argv[1::2], sys.argv[2::2]): row(nm, read_mrc(p_), (p_.replace('.mrc', '_even.mrc'), p_.replace('.mrc', '_odd.mrc')))
if len(sys.argv) > 1 and '--sharpen' in sys.argv:
    k = sys.argv.index('--sharpen'); src, dst = sys.argv[k + 1], sys.argv[k + 2]; v = read_mrc(src); F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); F *= np.exp(150. * s_inv ** 2 / 4.)
    try: res = fsc_res(read_mrc(src.replace('.mrc', '_even.mrc')), read_mrc(src.replace('.mrc', '_odd.mrc')))
    except Exception: res = 7.
    fc = 1. / res; F *= 0.5 * (1 + np.cos(np.pi * np.clip((s_inv - 0.85 * fc) / (0.3 * fc), 0, 1))); write_mrc(dst, np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F)))).astype(np.float32), SMPD); print('sharpened -> %s (low-pass %.1f A)' % (dst, res))
PYEOF
log "(1) ranking the $(grep -c . $FD/board/klem_latent.txt) board particles by the responsibility-weighted focal score"
$PY - <<PYEOF
import numpy as np
F = np.loadtxt('$FD/board/focal_scores.txt', comments='#'); fs = F[:, 1]; R = np.load('$FD/board/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); ngrid = len(R['i1'])
s = r[:, :ngrid] @ fs[:ngrid]; order = np.argsort(-s); N = len(np.load('$META')['flags']); st = np.zeros(N, int); N1 = $N1; N2 = $N2
st[pinds[order[:N1]] - 1] = 1; st[pinds[order[N1:N1 + N2]] - 1] = 2; st[pinds[order[-N1:]] - 1] = 3
np.savetxt('$D/states.txt', st, fmt='%d'); np.savetxt('$D/focal_score_per_particle.txt', np.c_[pinds, s], fmt='%d %.4f', header='row focal_score')
q = np.percentile(s, [10, 25, 50, 75, 90, 95, 99]); print('per-particle focal score: quantiles 10/25/50/75/90/95/99 = ' + ' '.join('%+.3f' % x for x in q))
for t in (0., 0.2, 0.4, 0.6, 0.8): print('  particles with score > %+.1f: %d (%.1f %% of the board)' % (t, int((s > t).sum()), 100. * np.mean(s > t)))
print('top %d >= %+.3f, next %d >= %+.3f, bottom %d <= %+.3f' % (N1, s[order[N1 - 1]], N2, s[order[N1 + N2 - 1]], N1, s[order[-N1]]))
sel1 = np.zeros(N, int); sel1[pinds[order[:N1]] - 1] = 1; np.savetxt('$D/states_top.txt', sel1, fmt='%d')
PYEOF
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=3 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for w in 01 02 03; do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps/ 2>/dev/null; done; done
$PY $SC $D/maps/recvol_state01.mrc "top $N1 by focal score (family poses)" $D/maps/recvol_state02.mrc "next $N2" $D/maps/recvol_state03.mrc "bottom $N1 (undocked)" $KF/closed_focal_closedfam/rank/maps/recvol_state01.mrc "closed-family docked 15k (reference)"
$PY $D/vis.py $D/maps/recvol_state01.mrc "top $N1 (family poses)" $D/maps/recvol_state02.mrc "next $N2 (family poses)" $D/maps/recvol_state03.mrc "bottom $N1 (undocked)" $KF/closed_focal_closedfam/rank/maps/recvol_state01.mrc "closed-family docked 15k (reference)"
log "(2) lever-arm refinement of the top $N1 as one state: refine3D refine=neigh athres=${ATHRES:-8} deg, full 232 A mask, 10 -> ${LPSTOP:-6} A, maxits ${MAXITS:-6}, nparts ${NP:-6}"
mkdir -p $D/refine && cd $D/refine && cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states_top.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '/state=1( |$)/' > $D/refine/ptcl3D_before.txt
simple_exec prg=refine3D projfile=p.simple nstates=1 refine=neigh athres=${ATHRES:-8} lpstart=10 lpstop=${LPSTOP:-6} maxits=${MAXITS:-6} mskdiam=232 pgrp=c1 nparts=${NP:-6} nthr=8 > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 'ERROR\|THROW' refine.log | cut -c1-200)"; tail -5 refine.log; exit 1; }
grep -a -m3 -i "refine\b.*neigh\|REFINE.*NEIGH\|>>> REFINE" refine.log | head -3
RF=$(ls -td $D/refine/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RF/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '/state=1( |$)/' > $D/refine/ptcl3D_after.txt
$PY - <<PYEOF
import numpy as np
def rd(fn, keys):
    return np.array([[float(dict(t.split('=', 1) for t in l.split() if '=' in t).get(k, 0)) for k in keys] for l in open(fn) if l.strip() and not l.startswith('#')])
def rotmat(ang, axis):
    r = np.radians(ang); c, s = np.cos(r), np.sin(r); R = np.eye(3)
    if axis == 3: R[0, 0], R[0, 1], R[1, 0], R[1, 1] = c, -s, s, c
    else:         R[0, 0], R[0, 2], R[2, 0], R[2, 2] = c, -s, s, c
    return R
def e2m(e): return rotmat(e[2], 3) @ rotmat(e[1], 2) @ rotmat(e[0], 3)
b = rd('$D/refine/ptcl3D_before.txt', ['e1', 'e2', 'e3', 'x', 'y']); a = rd('$D/refine/ptcl3D_after.txt', ['e1', 'e2', 'e3', 'x', 'y']); n = min(len(a), len(b))
ang = np.array([np.degrees(np.arccos(np.clip((np.trace(e2m(a[i, :3]) @ e2m(b[i, :3]).T) - 1) / 2, -1, 1))) for i in range(n)]); sh = np.hypot(a[:n, 3] - b[:n, 3], a[:n, 4] - b[:n, 4])
print('pose change by the lever-arm refinement over %d particles: angle median %.2f deg (p90 %.2f, max %.1f); shift median %.2f px (p90 %.2f)' % (n, np.median(ang), np.percentile(ang, 90), ang.max(), np.median(sh), np.percentile(sh, 90)))
PYEOF
log "(3) unmasked reconstruction at the lever-arm poses"; mkdir -p $D/refine/rec && cd $D/refine/rec && cp -f $RF/p.simple p.simple
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p maps; for h in "" _even _odd; do cp -f $RR/recvol_state01$h.mrc maps/ 2>/dev/null; done; cd $D
$PY $SC $D/maps/recvol_state01.mrc "top $N1, family poses (before)" $D/refine/rec/maps/recvol_state01.mrc "top $N1, lever-arm poses (after)"
mkdir -p $D/for_download; $PY $D/vis.py $D/maps/recvol_state01.mrc "top $N1, family poses (before)" $D/refine/rec/maps/recvol_state01.mrc "top $N1, lever-arm poses (after)" --sharpen $D/refine/rec/maps/recvol_state01.mrc $D/for_download/conf18_docked_lever_sharpB150.mrc
$PY $D/vis.py --sharpen $D/maps/recvol_state01.mrc $D/for_download/conf18_docked_top${N1}_sharpB150.mrc | tail -1
log DONE
