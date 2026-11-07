#!/bin/bash
#SBATCH -J rnasep_region_align_rec
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Node policy, second half, on the purest docked set (the within-10k): klem's strictly local pose step against the
# cell reference masked to the P12/P10.1 REGION (and, as a control, to a generous CORE = envelope minus the region) on
# a one-cell board of these particles at the family poses; the moved poses imported into a fresh project (stack order
# via stkind/indstk, CTF from the meta) and reconstructed unmasked by reconstruct3D (gold-standard halves), scored for
# the contact against the family-pose map of the same particles. env: SUBSET, STEP (3), PRIOR (4), EVERY (2), CTFFLAG (flip: the denoised stacks are phase-flipped)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=$A/rnasep_klem/expand1; P=$E/fam1_pipeline; KD=$P/klem_fb; META=$P/meta_fb.npz
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt; STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
FD=$KD/closed_focal; RK=$KD/closed_focal_docked30k/rank; SUBSET=${SUBSET:-$RK/subset_top10k.txt}; D=$KD/region_align_${TAG:-10k}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') ra] $*"; }; VIS=$FD/lever/vis.py
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation
reg = read_mrc('$FD/focal_mask_p12.mrc'); cons = read_mrc('$CONS'); g = gaussian_filter(cons, 2.0); dens = g > g.mean() + 0.5 * g.std()
core = dens & ~binary_dilation(reg > 0.5, iterations=3); soft = gaussian_filter(core.astype(float), 2.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1)
write_mrc('$D/core2_mask.mrc', soft.astype(np.float32), 1.464); print('masks: generous envelope %d voxels, core2 (envelope minus the dilated region) %d, region %d' % (dens.sum(), (soft > 0.5).sum(), (reg > 0.5).sum()))
PYEOF
export OMP_NUM_THREADS=1
run(){ tag=$1; shift; log "klem one-cell board of $(grep -c . $SUBSET) particles, pose step: $tag"; $PY $W/klem.py $STK none none $D/$tag --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 1 --n2 1 --nfree 0 --nproc 16 --init random --nsub 0 --subset $SUBSET --ctf auto \
    --lps "10 10 10 8 8 8" --sigmas "0.7 0.7 0.7 0.7 0.7 0.7" --mskdiam 232 --tag "ra $tag" "$@" 2>&1 | grep "subset\|pose step (\|Error\|Traceback" | cut -c1-200
  $PY $W/closed_arm_cells.py $D/$tag --top 1 2>&1 | grep -v Warning | sed "s/^/[$tag] /"; return 0; }
[ -s $D/lever/ptcl3D_pose_it6.txt ] || run lever --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell
[ -s $D/region/ptcl3D_pose_it6.txt ] || run region --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell --pose_mask $FD/focal_mask_p12.mrc
[ -s $D/core2/ptcl3D_pose_it6.txt ] || run core2  --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell --pose_mask $D/core2_mask.mrc
export OMP_NUM_THREADS=8; ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort > $D/stktab.txt
for tag in unmoved lever region core2; do
  [ -s $D/$tag/maps/recvol_state01.mrc ] && { log "$tag: map exists, skipped"; continue; }
  mkdir -p $D/$tag; $PY - <<PYEOF
import numpy as np, glob
M = np.load('$META'); N = len(M['flags']); e1, e2, e3, x, y = (M[k].astype(float).copy() for k in ('e1', 'e2', 'e3', 'x', 'y')); dfx, dfy, ang = (M[k].astype(float) for k in ('dfx', 'dfy', 'angast')); stk = M['stkind'].astype(int); ind = M['indstk'].astype(int); eo = M['eo'].astype(int) if 'eo' in M else np.arange(N) % 2
fs = sorted(glob.glob('$D/$tag/ptcl3D_pose_it*.txt'), key=lambda f: int(f.split('_it')[-1].split('.')[0])); rows = [dict(t.split('=', 1) for t in l.split()) for l in open(fs[-1]) if l.strip()] if fs else []
pind = np.array([int(r['pind']) for r in rows]) - 1 if rows else np.loadtxt('$SUBSET', dtype=int) - 1; state = np.zeros(N, int); state[pind] = 1
d0 = np.array([np.degrees(np.arccos(np.clip(np.cos(np.radians(float(r['e2']) - e2[p])), -1, 1))) for r, p in zip(rows[:2000], pind[:2000])]) if rows else np.zeros(1)
for r, p in zip(rows, pind): e1[p], e2[p], e3[p] = float(r['e1']), float(r['e2']), float(r['e3'])
order = np.lexsort((ind, stk)); n_stk = [int(np.fromfile(f.strip(), np.int32, 256)[2]) for f in open('$D/stktab.txt')]
ok = all(np.array_equal(ind[order][stk[order] == k], np.arange(1, n_stk[k - 1] + 1)) for k in range(1, len(n_stk) + 1)) and len(order) == sum(n_stk)
print('%s: %d rows selected, %d moved (%s; |d e2| median vs the family poses %.2f deg); table in stack order: %d rows, contiguous indstk per stack: %s' % ('$tag', int(state.sum()), len(rows), fs[-1].split('/')[-1] if fs else 'unmoved', np.median(d0), len(order), ok))
with open('$D/$tag/oritab.txt', 'w') as f:
    for i in order: f.write('e1=%.4f e2=%.4f e3=%.4f x=%.3f y=%.3f dfx=%.5f dfy=%.5f angast=%.3f eo=%d state=%d\n' % (e1[i], e2[i], e3[i], x[i], y[i], dfx[i], dfy[i], ang[i], eo[i], state[i]))
np.savetxt('$D/$tag/states_stackorder.txt', state[order], fmt='%d')
PYEOF
  cd $D/$tag && rm -rf proj && simple_exec prg=new_project projname=proj > np.log 2>&1 && cd proj || { log "new_project FAILED ($tag)"; continue; }
  simple_exec prg=import_particles projfile=proj.simple stktab=$D/stktab.txt oritab=$D/$tag/oritab.txt ctf=${CTFFLAG:-flip} smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED ($tag): $(grep -a -m2 'ERROR' import.log | cut -c1-200)"; tail -3 import.log; cd $D; continue; }
  simple_exec prg=selection projfile=proj.simple oritype=ptcl3D infile=$D/$tag/states_stackorder.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED ($tag)"; cd $D; continue; }
  log "import + selection done ($tag): state-1 rows $(simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null | grep -c ' state=1\.') ; eo values: $(simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null | grep -o ' eo=[^ ]*' | sort | uniq -c | tr '\n' ' ')"
  simple_exec prg=reconstruct3D projfile=proj.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED ($tag): $(grep -a -m2 ERROR rec.log | cut -c1-200)"; cd $D; continue; }
  RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $D/$tag/maps; for h in "" _even _odd; do cp -f $RR/recvol_state01$h.mrc $D/$tag/maps/ 2>/dev/null; done; cd $D
done
$PY $VIS ${REFMAP:-$RK/maps/recvol_state01.mrc} "${TAG:-10k}, family poses" $D/unmoved/maps/recvol_state01.mrc "re-imported unmoved (control)" $D/lever/maps/recvol_state01.mrc "whole-reference-aligned poses (lever)" $D/region/maps/recvol_state01.mrc "REGION-aligned poses" $D/core2/maps/recvol_state01.mrc "core-aligned poses"
log DONE
