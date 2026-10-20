#!/bin/bash
#SBATCH -J toy10w_fboard_pose
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The Fortran board's whole-reference pose step on the wide toy with the truth known: one-cell board of the TRUE
# conformer-18 particles at the consensus-refined poses (error 3.99 deg median), pose_every 2 / step 3 / prior 4
# over 6 rounds (the Python test gave 3.99 -> 2.89 deg, arm 0.80 -> 0.83); pose error vs truth from klem_poses.txt,
# closed-arm occupancy of the cell map. env: BDIR, DONEFILE
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${BDIR:-$HOME/SIMPLE_klem2/build}; DONE=${DONEFILE:-make9.done}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; S=$W/simw; L=$W/latw; T=$L/fboard_pose_c18; mkdir -p $T; log(){ echo "[$(date '+%F %T') fp] $*"; }
n=0; while [ $n -lt 90 ] && [ ! -f $B/$DONE ]; do sleep 60; n=$((n+1)); done; cat $B/$DONE; grep -q "_EXIT=0" $B/$DONE || { log "build failed"; exit 1; }
Z=$W/lat/klem_dump20k/klem_estep_ref.npz
$PY - <<PYEOF
import numpy as np
Z = np.load('$Z'); sub = np.loadtxt('$L/conf18_subset.txt', dtype=int)
for tag, extra in (('baseline', 'pose_every 0'), ('lever', 'pose_every 2\npose_step 3\npose_prior 4')):
    open('$T/board_%s.txt' % tag, 'w').write('n1 1\nn2 1\nnfree 0\nlps 10 10 10 8 8 8\nsigmas 0.7 0.7 0.7 0.7 0.7 0.7\nshrink 0.02\neps0 0.05\nbalance 0\nfocus_mask none\nsubset $L/conf18_subset.txt\nnsub 0\nseed 1\nss %g\nfs %g\nconv %d\nsign %g\nctfabs 0\nwrite_halves 0\n%s\n' % (float(Z['ss']), float(Z['fs']), int(Z['conv']), float(Z['sign']), extra))
print('configs written for %d true conformer-18 particles' % len(sub))
PYEOF
PJ=$(ls -td $W/projw/*_refine3D | head -1)/projw.simple
for tag in baseline lever; do
  log "klem_board $tag"; t0=$(date +%s)
  simple_exec prg=klem_board projfile=$PJ infile=$T/board_$tag.txt dir=$T/$tag nthr=16 mkdir=no > $T/$tag.log 2>&1 || { log "klem_board FAILED ($tag): $(grep -a -m3 'ERROR\|STOP\|Error' $T/$tag.log | cut -c1-200)"; tail -5 $T/$tag.log; continue; }
  grep -a ">>> it\|>>> pose\|>>> loading" $T/$tag.log | cut -c1-180; log "$tag done in $(( $(date +%s) - t0 )) s"
  $PY $W/closed_arm_cells.py $T/$tag --top 1 2>&1 | grep -v Warning | sed "s/^/[$tag] /"
  [ -s $T/$tag/klem_poses.txt ] && $PY - <<PYEOF
import numpy as np
def rd(fn, keys): return np.array([[float(dict(t.split('=', 1) for t in l.split() if '=' in t).get(k, 0)) for k in keys] for l in open(fn) if l.strip() and not l.startswith('#')])
def rotmat(ang, axis):
    r = np.radians(ang); c, s = np.cos(r), np.sin(r); R = np.eye(3)
    if axis == 3: R[0, 0], R[0, 1], R[1, 0], R[1, 1] = c, -s, s, c
    else: R[0, 0], R[0, 2], R[2, 0], R[2, 2] = c, -s, s, c
    return R
def e2m(e): return rotmat(e[2], 3) @ rotmat(e[1], 2) @ rotmat(e[0], 3)
t = rd('$S/all.txt', ['e1', 'e2', 'e3']); P = np.loadtxt('$T/$tag/klem_poses.txt', comments='#'); r0 = rd('$L/ptcl3D_r0.txt', ['e1', 'e2', 'e3'])
ang = []; ang0 = []
for row in P:
    i = int(row[0]) - 1; ang.append(np.degrees(np.arccos(np.clip((np.trace(e2m(row[1:4]) @ e2m(t[i]).T) - 1) / 2, -1, 1)))); ang0.append(np.degrees(np.arccos(np.clip((np.trace(e2m(r0[i]) @ e2m(t[i]).T) - 1) / 2, -1, 1))))
print('[$tag] pose error vs truth: start median %.2f deg (p90 %.2f) -> after the Fortran pose step median %.2f deg (p90 %.2f) on %d particles' % (np.median(ang0), np.percentile(ang0, 90), np.median(ang), np.percentile(ang, 90), len(ang)))
PYEOF
done
log DONE
