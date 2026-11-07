#!/bin/bash
#SBATCH -J rnasep_import_check2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Which stacks does the family-pose project reconstruct from? Its stk segment listed; then the within-10k set re-imported
# with THOSE stacks (unmoved poses) and reconstructed: must reproduce the family-pose map (contact 1.73) if the stacks
# were the only difference (the denoised stacks are what klem classifies on; the maps may need the raw particles).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; P=$A/rnasep_klem/expand1/fam1_pipeline; KD=$P/klem_fb; META=$P/meta_fb.npz; RK=$KD/closed_focal_docked30k/rank
D=$KD/region_align_10k/import_check2; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') ic2] $*"; }; VIS=$KD/closed_focal/lever/vis.py; SUBSET=$RK/subset_top10k.txt
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple
simple_exec prg=print_project_field projfile=$PS oritype=stk mkdir=no 2>/dev/null > $D/family_stk.txt; echo "== family project stk segment ($(grep -c . $D/family_stk.txt) entries)"; grep -o "stk=/[^ ]*" $D/family_stk.txt | head -3; grep -o "ctf=[^ ]*\|ctfflag=[^ ]*\|fraca=[^ ]*\|kv=[^ ]*\|cs=[^ ]*\|phshift=[^ ]*" $D/family_stk.txt | sort | uniq -c | head -12
grep -o "stk=/[^ ]*" $D/family_stk.txt | sed 's/stk=//' > $D/stktab_family.txt; CTFFLAG=$(grep -o "ctf=[^ ]*" $D/family_stk.txt | head -1 | sed 's/ctf=//'); FRACA=$(grep -o "fraca=[^ ]*" $D/family_stk.txt | head -1 | sed 's/fraca=//'); echo "import with ctf=${CTFFLAG:-yes} fraca=${FRACA:-0.1}"
echo "== denoised stktab for comparison"; ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | head -3
$PY - <<PYEOF
import numpy as np
M = np.load('$META'); N = len(M['flags']); e1, e2, e3, x, y = (M[k].astype(float) for k in ('e1', 'e2', 'e3', 'x', 'y')); dfx, dfy, ang = (M[k].astype(float) for k in ('dfx', 'dfy', 'angast')); stk = M['stkind'].astype(int); ind = M['indstk'].astype(int)
state = np.zeros(N, int); state[np.loadtxt('$SUBSET', dtype=int) - 1] = 1; order = np.lexsort((ind, stk))
with open('$D/oritab.txt', 'w') as f:
    for i in order: f.write('e1=%.4f e2=%.4f e3=%.4f x=%.3f y=%.3f dfx=%.5f dfy=%.5f angast=%.3f state=%d\n' % (e1[i], e2[i], e3[i], x[i], y[i], dfx[i], dfy[i], ang[i], state[i]))
print('oritab: %d rows, %d state 1' % (N, state.sum()))
PYEOF
rm -rf proj && simple_exec prg=new_project projname=proj > np.log 2>&1 && cd proj || { log "new_project FAILED"; exit 1; }
simple_exec prg=import_particles projfile=proj.simple stktab=$D/stktab_family.txt oritab=$D/oritab.txt ctf=${CTFFLAG:-yes} smpd=1.464 kv=300 cs=2.7 fraca=${FRACA:-0.1} mkdir=no > import.log 2>&1 || { log "import FAILED: $(grep -a -m2 ERROR import.log | cut -c1-200)"; tail -3 import.log; exit 1; }
simple_exec prg=reconstruct3D projfile=proj.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $D/maps; for h in "" _even _odd; do cp -f $RR/recvol_state01$h.mrc $D/maps/ 2>/dev/null; done; cd $D
$PY $VIS $RK/maps/recvol_state01.mrc "within-10k, family-pose project" $D/maps/recvol_state01.mrc "re-imported with the project CTF flag, unmoved poses" $KD/region_align_10k/import_check/plain/maps/recvol_state01.mrc "re-imported with the denoised stacks (before)"
log DONE
