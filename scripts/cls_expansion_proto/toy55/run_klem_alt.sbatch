#!/bin/bash
#SBATCH -J toy10_klem_alt
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Alternation rounds on toy10: klem (50k EM + final E-step over all 200k) -> 7 bins -> refine3D_states (local, 20->8 A,
# starting from the previous round's poses) -> next round. env: NROUND (2), KDIR (start latent dir, default lat/klem_junk),
# PROJSRC (project with the starting poses; default the consensus refine3D project), ALTDIR (output root)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; NB=${NB:-7}; NP=${NP:-10}; T=$W/sim; OUT=$W/lat; NROUND=${NROUND:-2}
RP=$(ls -td $W/proj/*_refine3D | head -1); PROJSRC=$(ls $RP/*.simple | head -1); POSETAB=${POSETAB:-}; KDIR=${KDIR:-$OUT/klem_junk}; ALT=${ALTDIR:-$OUT/alt}; mkdir -p $ALT
LOG=$(ls -t $W/toy10_battery.*.out | head -1); AXIS=$(grep -o "truth axis (PC1 of the conformers, STATES order): .*" $LOG | head -1 | sed 's/.*: //'); STATES=1,10,18,27,28,37,46,48,49,51; FD=$(ls -td $OUT/flex12/*_flex_pca | head -1)
log(){ echo "[$(date '+%F %T') alt] $*"; }
for R in $(seq 1 $NROUND); do D=$ALT/round$R; mkdir -p $D && cd $D; export OMP_NUM_THREADS=8
  if [ "${STATEMODE:-bins}" = cells ]; then   # the model's own states (argmax cell of every on-grid particle) as the refinement references
    $PY - <<PYEOF
import numpy as np
L = np.loadtxt('$KDIR/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); gm = L[:, 5]; ok = (gm > 0.5) & (am <= $NB)
st = np.zeros(199996, int); st[pind[ok] - 1] = am[ok]; np.savetxt('$D/states.txt', st, fmt='%d'); print('states = argmax cells: %d particles over %d states (min %d, max %d per state)' % (ok.sum(), len(np.unique(am[ok])), np.bincount(am[ok])[1:].min(), np.bincount(am[ok]).max()))
PYEOF
  else $PY $W/klem_bins.py $KDIR none $D/states.txt --nb $NB --nptcl 199996; fi
  if [ -n "$POSETAB" ]; then   # fresh project from the previous round's refined poses (a copied refine3D_states project loses its relative paths)
    mkdir -p $D/proj_in && cd $D/proj_in && simple_exec prg=new_project projname=p > np.log 2>&1 && cd p && simple_exec prg=import_particles projfile=p.simple stk=$T/all.mrcs oritab=$POSETAB ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 \
      && simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=48000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1 && PROJSRC=$D/proj_in/p/p.simple || { log "re-import FAILED"; exit 1; }
    cd $D; fi
  cp -f $PROJSRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1
  simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/seed_ptcl3D.txt; log "round $R seed: $($PY $W/pose_error.py $D/seed_ptcl3D.txt $T/all.txt --tag poses)"
  simple_exec prg=refine3D_states projfile=p.simple nstates=$NB pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-6} mskdiam=200 pgrp=c1 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
  RD=$(ls -td $D/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $D/refined_ptcl3D.txt
  log "round $R after refine3D_states: $($PY $W/pose_error.py $D/refined_ptcl3D.txt $T/all.txt --tag poses)"
  awk '{ gsub(/state=[0-9.]+/, "state=1"); print }' $D/refined_ptcl3D.txt > $D/refined_state1.txt; export OMP_NUM_THREADS=1; mkdir -p $D/klem
  $PY $W/klem.py $T/all.mrcs $D/refined_state1.txt $T/gt_particles.txt $D/klem --cons $RP/recvol_state01.mrc --flexdir "$FD" --states $STATES --axis=$AXIS --n1 7 --n2 7 --nfree 2 --nproc 16 --init pca --nsub 50000 --final_all --tag "alt round $R" 2>&1 | grep --line-buffered -v "^   [0-9]* / " | grep -E "final|it 6:|it 9:"
  KDIR=$D/klem; POSETAB=$D/refined_state1.txt
done
log DONE
