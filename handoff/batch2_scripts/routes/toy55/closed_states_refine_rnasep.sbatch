#!/bin/bash
#SBATCH -J rnasep_closed_states_refine
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The suspected-closed states REFINED AS STATES (user: not a resolution matter): the closed family (66k) with three
# states from the focal ranking (1 = docked top 15k, 2 = next 10k, 3 = the rest) refined by refine3D_states with the
# LOCAL pose policy (the mode that produced the family poses; plain refine3D ignores pose_policy), 10 -> LPSTOP A,
# wide mask; every state then reconstructed unmasked and scored for the contact; the docked map sharpened gently.
# env: LPSTOP (6), MAXITS (8), NP (10), RTAG
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; RK=$KF/closed_focal_closedfam/rank; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
D=$KF/closed_states_refine${RTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') csr] $*"; }
$PY - <<PYEOF
import numpy as np
st = np.loadtxt('$RK/states.txt', dtype=int); sub = np.loadtxt('$KF/closed_zoom/closed_subset.txt', dtype=int); N = len(st); lab = np.zeros(N, int); lab[sub - 1] = 3; lab[st == 1] = 1; lab[st == 2] = 2
np.savetxt('$D/states.txt', lab, fmt='%d'); print('closed family: %d particles -> states docked %d, intermediate %d, rest %d' % ((lab > 0).sum(), (lab == 1).sum(), (lab == 2).sum(), (lab == 3).sum()))
PYEOF
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
log "refine3D_states nstates=3 pose_policy=local 10 -> ${LPSTOP:-6} A maxits ${MAXITS:-8} mskdiam 232"
simple_exec prg=refine3D_states projfile=p.simple nstates=3 pose_policy=local lpstart=10 lpstop=${LPSTOP:-6} maxits=${MAXITS:-8} mskdiam=232 pgrp=c1 nparts=${NP:-10} nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; tail -5 refine.log; exit 1; }
RD=$(ls -td $D/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/refined_ptcl3D.txt
$PY - <<PYEOF
import numpy as np
st = []
for l in open('$D/refined_ptcl3D.txt'):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t); st.append(int(float(kv.get('state', 0))))
st = np.array(st); s0 = np.loadtxt('$D/states.txt', dtype=int); print('after refinement: state sizes %s (seed %s); particles that changed state %.2f' % (np.bincount(st, minlength=4)[1:], np.bincount(s0, minlength=4)[1:], float(np.mean(st[s0 > 0] != s0[s0 > 0]))))
PYEOF
log "unmasked reconstruction of the three refined states"; mkdir -p rec && cd rec && cp -f $RD/p.simple p.simple
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=3 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p maps; for w in 01 02 03; do for h in "" _even _odd; do cp -f $RR/recvol_state$w$h.mrc maps/ 2>/dev/null; done; done; cd $D
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464; n = 160; S = Shells(n, SMPD); sh = S.sh; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); rA = rr * SMPD; s_inv = rr / (n * SMPD)
lp = lambda v, res: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))) * (rr <= n * SMPD / res)))))
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); c27 = read_mrc(glob.glob(W + '/maps160/docked/conf27_*.mrc')[0])
g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); tip = arm & (rA > 95); ref = L18[arm].mean() / L18[core].mean()
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
print('contact VISIBILITY at 8 A (core median = 1) and arm occupancy:  %-30s %8s %8s %9s %6s %6s' % ('map', 'tip p99', 'tip p90', 'tip>0.5', 'arm', 'FSC'))
for name, vv, ev in (('deposited 18 (closed)', c18, None), ('deposited 27 (open)', c27, None), ('docked 15k, plain (before)', read_mrc('$RK/maps/recvol_state01.mrc'), ('$RK/maps/recvol_state01_even.mrc', '$RK/maps/recvol_state01_odd.mrc')), ('state 1 docked, REFINED as state', read_mrc('$D/rec/maps/recvol_state01.mrc'), ('$D/rec/maps/recvol_state01_even.mrc', '$D/rec/maps/recvol_state01_odd.mrc')), ('state 2 intermediate, refined', read_mrc('$D/rec/maps/recvol_state02.mrc'), ('$D/rec/maps/recvol_state02_even.mrc', '$D/rec/maps/recvol_state02_odd.mrc')), ('state 3 rest, refined', read_mrc('$D/rec/maps/recvol_state03.mrc'), ('$D/rec/maps/recvol_state03_even.mrc', '$D/rec/maps/recvol_state03_odd.mrc'))):
    L = lp(vv, 8.); L = L / max(np.median(L[core]), 1e-30); Lv = g(vv); occ = (Lv[arm].mean() / Lv[core].mean()) / ref
    try: res = fsc_res(read_mrc(ev[0]), read_mrc(ev[1])) if ev else float('nan')
    except Exception: res = float('nan')
    print('  %-54s %8.2f %8.2f %9.2f %6.2f %6.1f' % (name, np.percentile(L[tip], 99), np.percentile(L[tip], 90), np.mean(L[tip] > 0.5), occ, res))
v = read_mrc('$D/rec/maps/recvol_state01.mrc'); F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); F *= np.exp(150. * s_inv ** 2 / 4.)
try: res = fsc_res(read_mrc('$D/rec/maps/recvol_state01_even.mrc'), read_mrc('$D/rec/maps/recvol_state01_odd.mrc'))
except Exception: res = 7.
fc = 1. / res; F *= 0.5 * (1 + np.cos(np.pi * np.clip((s_inv - 0.85 * fc) / (0.3 * fc), 0, 1))); write_mrc('$D/conf18_docked_refined_state_sharpB150.mrc', np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F)))).astype(np.float32), SMPD); print('sharpened docked state map written (low-pass %.1f A)' % res)
PYEOF
log DONE
