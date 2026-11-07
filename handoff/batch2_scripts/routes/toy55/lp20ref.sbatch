#!/bin/bash
#SBATCH -J lp20ref
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# States assigned at <= 20 A at one s30split2 node: the node's focused K-state refinement (rflex6split step 4: n.simple, flex_states.txt
# selection, the flex-state refs, its focus mask, prob_neigh geom maxits 8 box_crop 88, focus3) repeated with ONLY lp = 20 (was 10) into
# rnasep_klem/lp20_s30split2_<node>; nothing in the split session's node dir is written. Then the per-state gate (nodetype_perstate2.py)
# at 17-10 A (2 shells past the lp-20 assignment), 13-10 A (also past the flex <= 15 A starting labels) and 9-7 A, with truth deviation
# shares and the overlap with the lp-10 states, and the lp-10 states' truth deviation shares at 9-7 A for comparison. env NODE (r3|r4|r5)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; GT=$T/simswing30/gt_particles.txt; TM=$T/simswing30/maps; MD=220
NODE=${NODE:?}; D=$K_/rflex_s30split2_$NODE/node_$NODE; O=$K_/lp20_s30split2_$NODE; mkdir -p $O && cd $O
log(){ echo "[$(date '+%F %T') lp20ref $NODE] $*"; }
[ -s $D/n.simple ] && [ -s $D/flex_states.txt ] && [ -s $D/mask.mrc ] && [ -s $D/members.txt ] && [ -s $D/refined_states.txt ] || { log "inputs missing in $D"; exit 1; }
KS=$(grep -c . $D/members.txt); VOLS=""; for j in $(seq 1 $KS); do VOLS="$VOLS vol$j=$D/refs/ref_state$(printf %02d $j).mrc"; done
t0=$(date +%s)
if [ -z "$(ls $O/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cp -f $D/n.simple r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$D/flex_states.txt mkdir=no > rsel.log 2>&1 || { log "state selection FAILED"; exit 1; }
  log "focused refine3D, $KS states from $D/refs, lp 20 (the node's own run: lp 10)"
  simple_exec prg=refine3D projfile=r.simple $VOLS nstates=$KS refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=20 maxits=8 center=no rec_backend=pcg pcg_mskfile=$D/mask.mrc pcg_focus=yes nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=4 nthr=16 > ref.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }
fi
RD=$(ls -td $O/*_refine3D | head -1); RP=$RD/r.simple; [ -s $RP ] || RP=$O/r.simple
log "refine3D $(( $(date +%s) - t0 )) s, $(grep -a -c '>>> ITERATION' ref.log) iterations, matching lp $(grep -a 'MATCHING  LOW-PASS' ref.log | awk '{print $NF}' | sort -n | uniq | tr '\n' ' ')"
simple_exec prg=print_project_field projfile=$RP oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > $O/refined_states.txt
$PY $T/scorecard.py $D/refined_states.txt $GT --codes swing:200:1 --tag "lp-10 states" 2>&1 | grep -v Warning
$PY $T/scorecard.py $O/refined_states.txt $GT --codes swing:200:1 --tag "lp-20 states" 2>&1 | grep -v Warning
$PY -u $T/nodetype_perstate2.py $RD $O/refined_states.txt --bands 17:10,13:10,9:7 --mskdiam $MD --gt $GT --tmaps $TM --cmp $D/refined_states.txt --tag "[$NODE lp20]" 2>&1
$PY -u $T/nodetype_perstate2.py $(ls -td $D/*_refine3D | head -1) $D/refined_states.txt --bands 9:7 --mskdiam $MD --gt $GT --tmaps $TM --tag "[$NODE lp10]" 2>&1
log "DONE ($(( $(date +%s) - t0 )) s)"
