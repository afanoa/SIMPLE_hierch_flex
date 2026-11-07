#!/bin/bash
#SBATCH -J klem_p2arm
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=230G
#SBATCH -t 10:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# E-step references from the tree (problem 2): refine3D one iteration at a time over all particles and the 34 rflex_s30a leaves, the
# references (vol1..vol34; lp set, so refine3D matches volN itself) being ARM's estimate of every leaf map from the current leaf
# half-maps (tree_wiener5.py --maps, fixed leaf tree post/tree): raw = the half-means (what refine3D matches without this), own = each
# leaf's own per-shell Wiener, tree2 = the tree posterior mean, two regions (borrowing along the tree). Per iteration the truth readout
# (p2_metrics.py, tree_wiener5 scores of every estimator); the refinement's new half-maps feed the next. Resumable. env ARM, NIT (2), LP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=24 OPENBLAS_NUM_THREADS=24
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_s30a; S=$T/simswing30
ARM=${ARM:?}; NIT=${NIT:-2}; P=$K_/p2estep/prep; W=$K_/p2estep/$ARM; TREE=$R/post/tree/merge_tree.txt; mkdir -p $W && cd $W; log(){ echo "[$(date '+%F %T') p2 $ARM] $*"; }
[ -e $P/prep.done ] || { log "prep missing ($P)"; exit 1; }
[ -s $W/lab_it1.txt ] || cp -f $P/labels.txt $W/lab_it1.txt; [ -s $W/list_it1.txt ] || cp -f $P/list.txt $W/list_it1.txt; [ -s $W/proj_it1.simple ] || cp -f $P/proj_it1.simple $W/proj_it1.simple
for it in $(seq 1 $NIT); do
  D=$W/it$it; mkdir -p $D; n=$((it + 1))
  [ -s $W/lab_it$n.txt ] && [ -s $W/list_it$n.txt ] && { log "iteration $it already done"; continue; }
  [ -s $D/tw/maps/$ARM/state_34.mrc ] || $PY $T/tree_wiener5.py $D/tw --list $W/list_it$it.txt --labels $W/lab_it$it.txt --tree $TREE --gt $S/gt_particles.txt --truth $S/maps --mskdiam 220 --maps $ARM 2>&1 | grep -v Warning | grep -v "smallest leaf" | sed "s/^/  it $it refs: /"
  [ -s $D/tw/maps/$ARM/state_34.mrc ] || { log "references FAILED it $it"; exit 1; }
  V=""; for k in $(seq 1 34); do V="$V vol$k=$D/tw/maps/$ARM/state_$(printf %02d $k).mrc"; done
  t0=$(date +%s)
  if [ -z "$(ls $D/*_refine3D/recvol_state34_odd_unfil.mrc 2>/dev/null)" ]; then cd $D; cp -f $W/proj_it$it.simple r.simple
    log "it $it: refine3D maxits=1, 34 states, lp ${LP:-10}, references $ARM"
    simple_exec prg=refine3D projfile=r.simple $V nstates=34 refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=220 lp=${LP:-10} maxits=1 center=no nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=2 nthr=24 > ref.log 2>&1 || { log "refine3D FAILED it $it: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }
    cd $W; fi
  RD=$(ls -td $D/*_refine3D | head -1); RP=$RD/r.simple; [ -s $RP ] || RP=$D/r.simple; rm -f $RD/*_part*_even.mrc $RD/*_part*_odd.mrc $RD/rho_*
  log "it $it: refine3D $(( ($(date +%s) - t0) / 60 )) min"
  simple_exec prg=print_project_field projfile=$RP oritype=ptcl3D mkdir=no 2>/dev/null > $D/dump.txt
  $PY $T/p2_metrics.py $D/dump.txt --gt $S/gt_particles.txt --true $S/all.txt --ref0 $W/lab_it1.txt --prev $W/lab_it$it.txt --out $D/labels.txt --tag "$ARM after it $it" 2>&1 | grep -v Warning
  awk -v d=$RD 'NR == FNR {c[$1]++; next} END {for (k = 1; k <= 34; k++) printf "%d %s/recvol_state%02d_even_unfil.mrc %s/recvol_state%02d_odd_unfil.mrc %d\n", k, d, k, d, k, c[k]}' $D/labels.txt /dev/null > $D/list.txt
  cp -f $RP $W/proj_it$n.simple && cp -f $D/list.txt $W/list_it$n.txt && cp -f $D/labels.txt $W/lab_it$n.txt
done
n=$((NIT + 1)); [ -s $W/final_tw/summary.txt ] || $PY $T/tree_wiener5.py $W/final_tw --list $W/list_it$n.txt --labels $W/lab_it$n.txt --tree $TREE --gt $S/gt_particles.txt --truth $S/maps --mskdiam 220 2>&1 | grep -v Warning | grep -v "smallest leaf" | sed "s/^/  final maps: /"
log DONE
