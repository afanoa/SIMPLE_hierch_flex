#!/bin/bash
#SBATCH -J toy10_kseed
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
log(){ echo "[$(date "+%F %T") kseed$TAG] $*"; }
FLEX=${FLEX:-$OUT/flex12}
[ -d $FLEX ] || { log "no flex dir $FLEX"; exit 1; }
AXIS=$(grep -o "truth axis (PC1 of the conformers, STATES order): .*" $(ls -t $W/toy10${TAG}_battery.*.out | head -1) | head -1 | sed 's/.*: //'); FD=$(ls -td $FLEX/*_flex_pca | head -1)
RP=$(ls -td $W/$PJ/*_refine3D | head -1); PROJFILE=$(ls $RP/*.simple | head -1); NPTCL=$(grep -vc "^#" $GT); LW="$PY $W/latent_windows.py"; SA="--states $STATES --axis=$AXIS"
D=$OUT/kseed${KTAG:-}; mkdir -p $D; log "START: flex $FD, K $K"; for g in 1 2 3 4; do FG=$(ls -td $OUT/flex_g$g/*_flex_pca 2>/dev/null | head -1); [ -n "$FG" ] && $PY $W/ptcl_latent_eval.py $FG $GT $SA --subset $OUT/members_g$g.txt --tag "4 local flex g$g" 2>&1 | grep -E "axis:|knn:|closed|Error|error" | sed "s/^/[toy10$TAG] /"; done
[ "${KSTART:-0}" = 0 ] && $PY $W/latent_kmeans_seed.py $FD $GT --nptcl $NPTCL --k $K --out $D/states_r1.txt $SA --ndim ${NDIM:-0} --tag "9 kseed${KTAG:-}" 2>&1 | sed "s/^/[toy10$TAG] /"
r=${KSTART:-0}; TABLE=$OUT/ptcl3D_r0.txt
for LP in ${KLPS:-20 14 10}; do r=$((r+1)); DR=$D/rec_r$r; mkdir -p $DR && cd $DR && cp -f $PROJFILE proj.simple
  simple_exec prg=selection projfile=proj.simple oritype=ptcl3D infile=$D/states_r$r.txt mkdir=no > sel.log 2>&1
  simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null > $DR/ptcl3D.txt
  simple_exec prg=reconstruct3D projfile=proj.simple pgrp=c1 mskdiam=200 nstates=$K nparts=$NP nthr=8 > rec.log 2>&1 || { log "reconstruct3D r$r FAILED"; break; }
  RD=$(ls -td $DR/*_reconstruct3D | head -1); mkdir -p $D/maps_r$r; for w in $(seq -f %02g 1 $K); do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps_r$r/ 2>/dev/null; done; done
  awk '!/^#/' $D/states_r$r.txt > /dev/null
  $LW estep $STK $DR/ptcl3D.txt $GT --mapdir $D/maps_r$r --nw $K --lp $LP --prev $D/states_r$r.txt --prior_sigma 1000 --out $D/estep_r$r.txt $SA --tag "9 kseed${KTAG:-} estep r$r lp$LP" 2>&1 | grep -v "^   [0-9]* / " | grep -v "smooth" | sed "s/^/[toy10$TAG] /"
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
log "DONE kseed$TAG -> $D"
