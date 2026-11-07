#!/bin/bash
#SBATCH -J rnasep_lever_pose_cont
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Lever-arm test done with a genuinely LOCAL refinement: refine3D refine=cont (continuous Levenberg-Marquardt
# pose refinement from the current pose, route ROUTE) on the docked top N1 as ONE state, full 232 A mask, 10 -> LPSTOP
# A. The refine=neigh run (53903807) re-searched a tenth of the particles by > 76 deg and shifts by 5.7 px median, so
# its arm verdict was void. Pose change (numeric state filter), unmasked reconstruction, arm/contact scores, sharpened map.
# env: FD, N1 (30000), ROUTE (shift_then_joint|joint), LPSTOP (6), MAXITS (6), NP (6), RTAG
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; FD=${FD:-$KF/closed_focal}
N1=${N1:-30000}; D=$FD/lever_pose_cont${RTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') lc] $*"; }; SC=$KF/closed_focal_closedfam/rank/score.py; VIS=$FD/lever/vis.py
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$FD/lever/states_top.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/ptcl3D_before_all.txt
log "refine3D refine=cont (snapshot has no pose_cont_route), one state (docked top $N1), full 232 A mask, 10 -> ${LPSTOP:-6} A, maxits ${MAXITS:-6}, nparts ${NP:-6}"
simple_exec prg=refine3D projfile=p.simple nstates=1 refine=cont lpstart=10 lpstop=${LPSTOP:-6} maxits=${MAXITS:-6} mskdiam=232 pgrp=c1 nparts=${NP:-6} nthr=8 > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 'ERROR\|THROW' refine.log | cut -c1-200)"; tail -5 refine.log; exit 1; }
grep -a -m4 "REFINEMENT MODE\|POSE_CONT\|pose_cont" refine.log | cut -c1-160
RF=$(ls -td $D/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RF/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/ptcl3D_after_all.txt
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
K = ['e1', 'e2', 'e3', 'x', 'y', 'state']; b = rd('$D/ptcl3D_before_all.txt', K); a = rd('$D/ptcl3D_after_all.txt', K); n = min(len(a), len(b)); m = (b[:n, 5] > 0.5) & (a[:n, 5] > 0.5)
ang = np.array([np.degrees(np.arccos(np.clip((np.trace(e2m(a[i, :3]) @ e2m(b[i, :3]).T) - 1) / 2, -1, 1))) for i in np.where(m)[0]]); sh = np.hypot(a[:n, 3] - b[:n, 3], a[:n, 4] - b[:n, 4])[m]
print('pose change by pose_cont over %d docked particles: angle median %.2f deg (p90 %.2f, max %.1f); shift median %.2f px (p90 %.2f); moved > 10 deg %.1f %%' % (m.sum(), np.median(ang), np.percentile(ang, 90), ang.max(), np.median(sh), np.percentile(sh, 90), 100. * np.mean(ang > 10)))
PYEOF
log "unmasked reconstruction at the pose_cont poses"; mkdir -p $D/rec && cd $D/rec && cp -f $RF/p.simple p.simple
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p maps; for h in "" _even _odd; do cp -f $RR/recvol_state01$h.mrc maps/ 2>/dev/null; done; cd $D
$PY $SC $FD/lever/maps/recvol_state01.mrc "top $N1, family poses (before)" $D/rec/maps/recvol_state01.mrc "top $N1, pose_cont poses (after)"
$PY $VIS $FD/lever/maps/recvol_state01.mrc "top $N1, family poses (before)" $D/rec/maps/recvol_state01.mrc "top $N1, pose_cont poses (after)"
mkdir -p $D/for_download; $PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc, Shells, fftc
SMPD = 1.464; n = 160; S = Shells(n, SMPD); sh = S.sh; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); s_inv = rr / (n * SMPD)
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
for src, dst in (('$D/rec/maps/recvol_state01.mrc', '$D/for_download/conf18_docked_pose_cont_sharpB150.mrc'),):
    v = read_mrc(src); F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); F *= np.exp(150. * s_inv ** 2 / 4.); res = fsc_res(read_mrc(src.replace('.mrc', '_even.mrc')), read_mrc(src.replace('.mrc', '_odd.mrc'))); fc = 1. / res
    F *= 0.5 * (1 + np.cos(np.pi * np.clip((s_inv - 0.85 * fc) / (0.3 * fc), 0, 1))); write_mrc(dst, np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F)))).astype(np.float32), SMPD); print('sharpened -> %s (low-pass %.1f A)' % (dst, res))
PYEOF
log DONE
