#!/bin/bash
#SBATCH -J cheapkids2
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Cheap view-balanced per-particle children of one rflex node (cheapkids.py: node flex latent centred per direction bin, PC1, K quantile
# bins inside every bin) -> reconstruct3D K states at the parent's poses (no refinement) -> truth: scorecard.py for the children and for
# the particle route's labels at the same node (flex states, refined states, coordinate leaves), subcls_mapscore.py for the children
# and the coordinate leaves. Times each stage. env RUN (rflex_s30a), NODE, K (10), TAG (swing30; '' = real data, no truth), MSKDIAM (220)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
RUN=${RUN:-rflex_s30a}; NODE=${NODE:?}; K=${K:-10}; TAG=${TAG-swing30}; MD=${MSKDIAM:-220}; GT=$T/sim$TAG/gt_particles.txt; TM=$T/sim$TAG/maps
RD=$K_/$RUN; ND=$RD/node_$NODE; O=$K_/cheapkids_${RUN#rflex_}_$NODE; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') cheapkids $RUN $NODE] $*"; }
set -- $(awk -v u=$NODE '$1 == u' $RD/queue.txt); PP=${3:-}; ROWS=${4:-}; FLEX=$(ls -td $ND/*_flex_pca | head -1)/flex_pca_coordinates.txt
[ -s "$PP" ] && [ -s "$ROWS" ] && [ -s "$FLEX" ] || { log "inputs missing"; exit 1; }
t0=$(date +%s)
[ -s $O/ptcl3D_parent.txt ] || simple_exec prg=print_project_field projfile=$PP oritype=ptcl3D mkdir=no > $O/ptcl3D_parent.txt 2> /dev/null
$PY $T/cheapkids.py --oris $O/ptcl3D_parent.txt --rows $ROWS --flex $FLEX --out $O --k $K $([ -n "$TAG" ] && echo --gt $GT) 2>&1 | grep -v Warning || { log "cheapkids.py FAILED"; exit 1; }
t1=$(date +%s); log "children: $((t1 - t0)) s (pose dump + bins + PC1)"
if [ -z "$(ls $O/rec/*_reconstruct3D/recvol_state$(printf %02d $K)_odd_unfil.mrc 2>/dev/null)" ]; then mkdir -p $O/rec && cd $O/rec
  cp -f $PP r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/children.txt mkdir=no > rsel.log 2>&1 \
    && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=6 mkdir=no > rup.log 2>&1 || { log "project FAILED"; exit 1; }
  simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=$K nparts=4 nthr=6 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; exit 1; }
  RR=$(ls -td $O/rec/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*; cd $O; fi
RR=$(ls -td $O/rec/*_reconstruct3D | head -1); t2=$(date +%s); log "reconstruct3D $K children: $((t2 - t1)) s"
[ -z "$TAG" ] && { log "DONE, no truth ($(( $(date +%s) - t0 )) s)"; exit 0; }
echo "== truth scorecards (same node members)"
$PY $T/scorecard.py $O/children.txt $GT --codes swing:200:1 --tag "cheap children (PC1, view-balanced)" 2>&1 | grep -v Warning
[ -s $ND/flex_states.txt ] && $PY $T/scorecard.py $ND/flex_states.txt $GT --codes swing:200:1 --tag "particle route: flex states" 2>&1 | grep -v Warning
[ -s $ND/refined_states.txt ] && $PY $T/scorecard.py $ND/refined_states.txt $GT --codes swing:200:1 --tag "particle route: focused-refined states" 2>&1 | grep -v Warning
[ -s $RD/coord_$NODE/leaf_states.txt ] && $PY $T/scorecard.py $RD/coord_$NODE/leaf_states.txt $GT --codes swing:200:1 --tag "particle route: coordinate leaves" 2>&1 | grep -v Warning
echo "== map scores"
$PY $T/subcls_mapscore.py $RR $O/children.txt $GT $TM --mskdiam $MD --tag cheap_children 2>&1 | grep -v Warning
[ -s $RD/coord_$NODE/leaf_states.txt ] && $PY $T/subcls_mapscore.py $RD/coord_$NODE/leaves/1_reconstruct3D $RD/coord_$NODE/leaf_states.txt $GT $TM --mskdiam $MD --tag coord_leaves 2>&1 | grep -v Warning
log "DONE ($(( $(date +%s) - t0 )) s)"
