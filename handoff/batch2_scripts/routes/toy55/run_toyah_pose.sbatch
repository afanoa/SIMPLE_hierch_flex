#!/bin/bash
#SBATCH -J toyah_pose
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Can the Fortran board's own local pose step replace refine3D_states at FAMILY scale (many states, 18 conformers +
# junk)? On the "ah" toy with the truth known: 7x7 board at the toy's (refined or perturbed) poses with the
# whole-reference pose step every 2 rounds (27-rotation grid, step 3 deg, prior 4 deg); pose error vs the truth before
# and after, over conformer particles, closed particles and per conformer. env: AHTAG (ah), EVERY (2), STEP (3), PRIOR (4), NT (16)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; AHTAG=${AHTAG:-ah}; S=$W/sim$AHTAG; L=$W/lat_$AHTAG; export TOY_SMPD=1.464; log(){ echo "[$(date '+%F %T') ahp] $*"; }
n=0; while [ $n -lt 600 ] && [ ! -s $L/ptcl3D_r0.txt ]; do sleep 60; n=$((n+1)); done; [ -s $L/ptcl3D_r0.txt ] || { log "no toy"; exit 1; }
RD=$(ls -td $W/proj$AHTAG/*_refine3D $W/proj$AHTAG/1_fastposes 2>/dev/null | head -1); PROJ=$RD/proj$AHTAG.simple; GT=$S/gt_particles.txt; NALL=$(grep -vc "^#" $GT); seq 1 $NALL > $L/subset_all.txt
export FB=${FB:-$HOME/SIMPLE_klem2/build} NTHR=${NT:-16}; O=$L/fP; log "board with the pose step (every ${EVERY:-2}, step ${STEP:-3}, prior ${PRIOR:-4}) on $NALL particles"
[ -s $O/klem_poses.txt ] || INIT=random NFREE=2 POSE_EVERY=${EVERY:-2} POSE_STEP=${STEP:-3} POSE_PRIOR=${PRIOR:-4} bash $W/fboard.sh $O $PROJ $L/subset_all.txt "14 14 12 12 10 10 8 8 8" "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" || exit 1
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
def angerr(P, T): return np.array([np.degrees(np.arccos(np.clip((np.trace(e2m(p) @ e2m(t).T) - 1) / 2, -1, 1))) for p, t in zip(P, T)])
t = rd('$S/all.txt', ['e1', 'e2', 'e3', 'state']); r0 = rd('$L/ptcl3D_r0.txt', ['e1', 'e2', 'e3']); n = min(len(t), len(r0)); conf = t[:n, 3].astype(int)
P = np.loadtxt('$O/klem_poses.txt', comments='#'); r1 = r0.copy(); pind = P[:, 0].astype(int) - 1; r1[pind[pind < n]] = P[pind < n, 1:4]
cl = set(int(v) for v in open('$S/closed_ids.txt').read().strip().split(',')); good = conf > 0; closed = np.array([c in cl for c in conf]); holo = conf >= 101
e0, e1 = angerr(r0[:n], t[:n, :3]), angerr(r1[:n], t[:n, :3])
print('POSE STEP at family scale (%d particles, %d conformer): pose error median %.2f -> %.2f deg (p90 %.2f -> %.2f); closed %.2f -> %.2f; holo %.2f -> %.2f; apo %.2f -> %.2f; junk %.1f -> %.1f' % (n, good.sum(), np.median(e0[good]), np.median(e1[good]), np.percentile(e0[good], 90), np.percentile(e1[good], 90), np.median(e0[closed]), np.median(e1[closed]), np.median(e0[holo]), np.median(e1[holo]), np.median(e0[good & ~holo]), np.median(e1[good & ~holo]), np.median(e0[~good]), np.median(e1[~good])))
print('  per conformer (before -> after, median deg): ' + ' '.join('%d:%.1f->%.1f' % (c, np.median(e0[conf == c]), np.median(e1[conf == c])) for c in np.unique(conf[good])))
print('  share of conformer particles improved %.2f, worsened by > 1 deg %.2f' % (np.mean(e1[good] < e0[good]), np.mean(e1[good] > e0[good] + 1)))
PYEOF
$PY $W/score_ah.py $O/klem_latent.txt $GT $S/closed_ids.txt --tag "board with pose step" 2>&1 | grep -v Warning | head -4
log DONE
