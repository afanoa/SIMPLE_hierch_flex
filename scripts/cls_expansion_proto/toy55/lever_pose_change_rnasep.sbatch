#!/bin/bash
#SBATCH -J rnasep_lever_pose_change
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Pose change of the lever-arm refinement (the in-job readout failed on the state=1.0000 filter): family poses vs the
# refine3D refine=neigh poses of the docked top N1, angle and shift statistics; both tables saved.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; FD=${FD:-$A/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal}; D=$FD/lever/refine; cd $D
RF=$(ls -td $D/*_refine3D | head -1); log(){ echo "[$(date '+%F %T') pc] $*"; }
simple_exec prg=print_project_field projfile=$D/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/ptcl3D_before_all.txt
simple_exec prg=print_project_field projfile=$RF/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/ptcl3D_after_all.txt
log "tables: before $(grep -c . $D/ptcl3D_before_all.txt) rows, after $(grep -c . $D/ptcl3D_after_all.txt) rows"
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
print('pose change by the lever-arm refinement over %d docked particles: angle median %.2f deg (p90 %.2f, max %.1f); shift median %.2f px (p90 %.2f); unchanged (< 0.1 deg) %.1f %%' % (m.sum(), np.median(ang), np.percentile(ang, 90), ang.max(), np.median(sh), np.percentile(sh, 90), 100. * np.mean(ang < 0.1)))
PYEOF
log DONE
