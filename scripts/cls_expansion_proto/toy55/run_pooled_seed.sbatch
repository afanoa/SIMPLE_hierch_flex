#!/bin/bash
#SBATCH -J toy10_pooled
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Step 7: states seeded from the 2D sub-class-POOLED flex latent (waits for the battery's pose-mean expansion), then
# hard EM against the particles: reconstruct3D (K states) -> cross-half E-step (flat prior) -> argmax labels, 3 rounds
# at 20, 14, 10 A, truth-scored every round. env: TAG (""), K (8), NP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; TAG=${TAG:-}; K=${K:-8}; NP=${NP:-10}; SIMDIR=${SIMDIR:-sim}; PJ=${PROJ:-proj}
T=$W/$SIMDIR; OUT=$W/lat$TAG; GT=$T/gt_particles.txt; STK=$T/all.mrcs; CSV=$W/${CSV:-rnasep_conformers10.csv}; STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $CSV)
log(){ echo "[$(date '+%F %T') pooled$TAG] $*"; }
n=0; while [ $n -lt 600 ] && [ ! -s $OUT/exppm/ptcl2D.txt ]; do sleep 60; n=$((n+1)); done
[ -s $OUT/exppm/ptcl2D.txt ] || { log "no pose-mean expansion table after 10 h"; exit 1; }
AXIS=$(grep -o "truth axis (PC1 of the conformers, STATES order): .*" $(ls -t $W/toy10${TAG}_battery.*.out | head -1) | head -1 | sed 's/.*: //'); FD=$(ls -td $OUT/flex/*_flex_pca | head -1)
RP=$(ls -td $W/$PJ/*_refine3D | head -1); PROJFILE=$(ls $RP/*.simple | head -1); NPTCL=$(grep -vc "^#" $GT); LW="$PY $W/latent_windows.py"; SA="--states $STATES --axis=$AXIS"
D=$OUT/pooled; mkdir -p $D; log "START: flex $FD, expansion $OUT/exppm/ptcl2D.txt, K $K"
$PY $W/pooled_seed.py $FD $OUT/exppm/ptcl2D.txt $GT --nptcl $NPTCL --k $K --out $D/states_r1.txt $SA --tag "7 pooled" 2>&1 | sed "s/^/[toy10$TAG] /"
r=0; TABLE=$OUT/ptcl3D_r0.txt
for LP in 20 14 10; do r=$((r+1)); DR=$D/rec_r$r; mkdir -p $DR && cd $DR && cp -f $PROJFILE proj.simple
  simple_exec prg=selection projfile=proj.simple oritype=ptcl3D infile=$D/states_r$r.txt mkdir=no > sel.log 2>&1
  simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null > $DR/ptcl3D.txt
  simple_exec prg=reconstruct3D projfile=proj.simple pgrp=c1 mskdiam=200 nstates=$K nparts=$NP nthr=8 > rec.log 2>&1 || { log "reconstruct3D r$r FAILED"; break; }
  RD=$(ls -td $DR/*_reconstruct3D | head -1); mkdir -p $D/maps_r$r; for w in $(seq -f %02g 1 $K); do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps_r$r/ 2>/dev/null; done; done
  awk '!/^#/' $D/states_r$r.txt > /dev/null
  $LW estep $STK $DR/ptcl3D.txt $GT --mapdir $D/maps_r$r --nw $K --lp $LP --prev $D/states_r$r.txt --prior_sigma 1000 --out $D/estep_r$r.txt $SA --tag "7 pooled estep r$r lp$LP" 2>&1 | grep -v "^   [0-9]* / " | grep -v "smooth" | sed "s/^/[toy10$TAG] /"
  [ -s $D/estep_r$r.txt ] || { log "E-step r$r failed"; break; }
  # next round's hard labels = argmax (column 3), particles outside keep 0
  $PY - $D/estep_r$r.txt $NPTCL $D/states_r$((r+1)).txt <<'PYEOF'
import sys, numpy as np
st = np.zeros(int(sys.argv[2]), int)
for l in open(sys.argv[1]):
    if l.startswith('#'): continue
    w = l.split(); st[int(w[0]) - 1] = int(w[2])
np.savetxt(sys.argv[3], st, fmt='%d'); print('labels -> %s (%d labelled, sizes %s)' % (sys.argv[3], (st > 0).sum(), np.bincount(st)[1:]))
PYEOF
done
log "DONE pooled$TAG -> $D"
