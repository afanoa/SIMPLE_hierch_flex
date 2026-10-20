#!/bin/bash
#SBATCH -J rnasep_import_check
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Why did the re-imported poses reconstruct wrongly (contact 0.30 for both variants)? (a) round trip of klem's pose
# dump convention on the meta angles; (b) the UNMOVED meta poses of the within-10k set written as an oritab in stack
# order, imported with stktab, reconstructed: must reproduce the family-pose map (contact 1.73); (c) the same with the
# shift signs flipped. Whichever variant reproduces the map is the right import convention.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KD=$P/klem_fb; META=$P/meta_fb.npz; RK=$KD/closed_focal_docked30k/rank
D=$KD/region_align_10k/import_check; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') ic] $*"; }; VIS=$KD/closed_focal/lever/vis.py; SUBSET=$RK/subset_top10k.txt
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); sys.path.insert(0, '$A/routes/pvbag/subclass_pool'); from view_tangent_toy import euler2m
def m2euler_simple(R):
    e2 = np.degrees(np.arccos(np.clip(R[2, 2], -1, 1)))
    if abs(R[2, 2]) < 1 - 1e-8: e1 = np.degrees(np.arctan2(-R[2, 1], R[2, 0])); e3 = np.degrees(np.arctan2(-R[1, 2], -R[0, 2]))
    else: e1 = np.degrees(np.arctan2(R[1, 0], R[0, 0])) * (1 if R[2, 2] > 0 else -1); e3 = 0.
    return np.array([e1 % 360., e2, e3 % 360.])
M = np.load('$META'); e1, e2, e3 = (M[k].astype(float) for k in ('e1', 'e2', 'e3')); rows = np.loadtxt('$SUBSET', dtype=int)[:3000] - 1; err = []
for i in rows:
    R = euler2m([-e1[i], e2[i], -e3[i]]); e = m2euler_simple(R); back = np.array([(-e[0]) % 360., e[1], (-e[2]) % 360.]); Rb = euler2m([-back[0], back[1], -back[2]])
    err.append(np.degrees(np.arccos(np.clip((np.trace(R @ Rb.T) - 1) / 2, -1, 1))))
print('(a) klem pose-dump round trip on %d meta rows: rotation error median %.4f deg, max %.4f' % (len(rows), np.median(err), np.max(err)))
PYEOF
ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort > $D/stktab.txt
for tag in plain flipshift; do
  $PY - <<PYEOF
import numpy as np
M = np.load('$META'); N = len(M['flags']); e1, e2, e3, x, y = (M[k].astype(float) for k in ('e1', 'e2', 'e3', 'x', 'y')); dfx, dfy, ang = (M[k].astype(float) for k in ('dfx', 'dfy', 'angast')); stk = M['stkind'].astype(int); ind = M['indstk'].astype(int)
sgn = -1. if '$tag' == 'flipshift' else 1.; state = np.zeros(N, int); state[np.loadtxt('$SUBSET', dtype=int) - 1] = 1; order = np.lexsort((ind, stk))
with open('$D/oritab_$tag.txt', 'w') as f:
    for i in order: f.write('e1=%.4f e2=%.4f e3=%.4f x=%.3f y=%.3f dfx=%.5f dfy=%.5f angast=%.3f state=%d\n' % (e1[i], e2[i], e3[i], sgn * x[i], sgn * y[i], dfx[i], dfy[i], ang[i], state[i]))
print('(%s) oritab written: %d rows, %d with state 1; sample row 0: e1 %.2f e2 %.2f e3 %.2f x %.2f y %.2f' % ('$tag', N, state.sum(), e1[order[0]], e2[order[0]], e3[order[0]], sgn * x[order[0]], sgn * y[order[0]]))
PYEOF
  mkdir -p $D/$tag && cd $D/$tag && rm -rf proj && simple_exec prg=new_project projname=proj > np.log 2>&1 && cd proj || { log "new_project FAILED"; cd $D; continue; }
  simple_exec prg=import_particles projfile=proj.simple stktab=$D/stktab.txt oritab=$D/oritab_$tag.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED ($tag): $(grep -a -m2 ERROR import.log | cut -c1-160)"; cd $D; continue; }
  simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk 'NR<=3' | cut -c1-200 | sed "s/^/[$tag imported row] /"
  simple_exec prg=reconstruct3D projfile=proj.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED ($tag)"; cd $D; continue; }
  RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $D/$tag/maps; for h in "" _even _odd; do cp -f $RR/recvol_state01$h.mrc $D/$tag/maps/ 2>/dev/null; done; cd $D
done
echo "== reference project row for comparison (family-pose project, first state-1 row)"; PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; simple_exec prg=print_project_field projfile=$PS oritype=ptcl3D mkdir=no 2>/dev/null | awk 'NR<=2' | cut -c1-220
$PY $VIS $RK/maps/recvol_state01.mrc "within-10k, family-pose project" $D/plain/maps/recvol_state01.mrc "re-imported, unmoved poses" $D/flipshift/maps/recvol_state01.mrc "re-imported, shifts flipped"
log DONE
