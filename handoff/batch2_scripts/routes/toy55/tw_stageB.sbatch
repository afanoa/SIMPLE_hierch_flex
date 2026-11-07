#!/bin/bash
#SBATCH -J klem_twstageB
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# tree-borrowed references inside a node's focused K-state refinement (proposal for rflex10 stage B; tested on a finished node): refine3D one
# iteration at a time (maxits=1) with the node refinement's own keys (prob_neigh geom, lp LP, box_crop 88, nspace 5000, pcg_focus with the node's
# mask), the references of every iteration = ARM's estimate of the states from the previous iteration's half-maps: raw = the state maps refine3D
# wrote (what its own next iteration would match), auto = the posterior mean over the node's states along their tree (state_merge.py --list,
# complete linkage -> tw_units.py --stage maps at the input box, no padding / low-pass): tree2 when the tree's off-diagonal R^2 >= 0.7, else emp2.
# Truth readout per iteration (p2_metrics.py) and of the final state maps (tree_wiener5.py). Start (START=flex): the node's flex-state labels and
# maps (rflex4 layout: n.simple, flex_states.txt, refs/, members.txt, mask.mrc); START=stageB: rflex10 layout (fstageB/b.simple, fstageB/*_reconstruct3D,
# rsmask.mrc). env RT, ND, ARMS ("raw auto"), NIT (4), LP (10), MSKDIAM (220), TAG (swing30: truth)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
RT=${RT:?}; ND=${ND:?}; N=$K_/rflex_$RT/node_$ND; NIT=${NIT:-4}; MD=${MSKDIAM:-220}; TAG=${TAG:-swing30}; S=$T/sim$TAG; O=$K_/tw_stageB/${RT}_$ND; mkdir -p $O
log(){ echo "[$(date '+%F %T') tw_stageB $RT $ND] $*"; }
if [ "${START:-flex}" = stageB ]; then P0=$N/fstageB/b.simple; MASK=$N/rsmask.mrc; RR0=$(ls -td $N/fstageB/*_reconstruct3D | head -1)
else MASK=$N/mask.mrc; fi
st(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null; }
labels(){ awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' $1; }
borrow(){ # $1 work dir, $2 list file (k even odd count) -> $1/refs.txt (vol keys) of the auto estimator
  mkdir -p $1 && cp -f $2 $1/units_list.txt
  [ -s $1/tree/merge_tree.txt ] || $PY $T/state_merge.py $1/tree --list $1/units_list.txt --link complete --ncut 2 > $1/tree.log 2>&1 || return 1
  local nb=$(head -c 4 $(awk 'NR == 1 {print $2}' $1/units_list.txt) | od -An -tu4 | tr -d ' ')
  $PY $T/tw_units.py $1 --stage maps --mapdir $1/unused --deliver none --box $nb --lp 0 --mskdiam $MD 2>&1 | grep -v Warning | grep "R^2\|FSC" | sed "s/^/   /"
  local est=$(awk '{print $1}' $1/choice.txt); awk -v d=$1/maps_$est '{printf " vol%d=%s/unit_%02d.mrc", $1, d, $1}' $1/units_list.txt > $1/refs.txt; echo "   references: $est"; }
for ARM in ${ARMS:-raw auto}; do W=$O/$ARM; mkdir -p $W
  # iteration-1 inputs: project with the start labels, start references (raw) or their borrowed version
  if [ "${START:-flex}" = stageB ]; then cp -f $P0 $W/proj_1.simple; KS=$(ls $RR0/recvol_state??.mrc | wc -l)
    awk -v d=$RR0 -v k=$KS 'BEGIN {for (s = 1; s <= k; s++) printf "%d %s/recvol_state%02d_even_unfil.mrc %s/recvol_state%02d_odd_unfil.mrc 0\n", s, d, s, d, s}' > $W/list_1.txt
    V1=$(awk -v d=$RR0 '{printf " vol%d=%s/recvol_state%02d.mrc", $1, d, $1}' $W/list_1.txt)
  else [ -s $W/proj_1.simple ] || { cp -f $N/n.simple $W/proj_1.simple && simple_exec prg=selection projfile=$W/proj_1.simple oritype=ptcl3D infile=$N/flex_states.txt mkdir=no > $W/sel.log 2>&1; }
    cp -f $N/members.txt $W/list_1.txt; V1=$(awk -v d=$N/refs '{printf " vol%d=%s/ref_state%02d.mrc", $1, d, $1}' $W/list_1.txt); fi
  st $W/proj_1.simple > $W/dump_0.txt; labels $W/dump_0.txt > $W/lab_0.txt
  [ "${START:-flex}" = stageB ] && awk 'NR == FNR {if ($1 > 0) c[$1]++; next} {print $1, $2, $3, c[$1] + 0}' $W/lab_0.txt $W/list_1.txt > $W/list_1.tmp && mv $W/list_1.tmp $W/list_1.txt
  [ -s $W/metrics_0.txt ] || $PY $T/p2_metrics.py $W/dump_0.txt --gt $S/gt_particles.txt --true $S/all.txt --ref0 $W/lab_0.txt --nsmall 3 --tag "$ARM start" 2>&1 | grep -v Warning > $W/metrics_0.txt; cat $W/metrics_0.txt
  for it in $(seq 1 $NIT); do D=$W/it$it; mkdir -p $D
    if [ -z "$(ls $D/*_refine3D/recvol_state01_odd_unfil.mrc 2>/dev/null)" ]; then
      if [ $ARM = raw ]; then [ $it = 1 ] && V=$V1 || V=$(cat $W/rawrefs_$it.txt)
      else borrow $D/borrow $W/list_$it.txt || { log "$ARM it $it borrowing FAILED"; exit 1; }; V=$(cat $D/borrow/refs.txt); fi
      KS=$(wc -l < $W/list_$it.txt); cd $D; cp -f $W/proj_$it.simple r.simple; t0=$(date +%s)
      simple_exec prg=refine3D projfile=r.simple $V nstates=$KS refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=${LP:-10} maxits=1 center=no rec_backend=pcg pcg_mskfile=$MASK pcg_focus=yes nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=4 nthr=12 > ref.log 2>&1 || { log "$ARM it $it refine3D FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }
      log "$ARM it $it: refine3D $(( ($(date +%s) - t0) / 60 )) min"; cd $O; fi
    RD=$(ls -td $D/*_refine3D | head -1); RP=$RD/r.simple; [ -s $RP ] || RP=$D/r.simple; rm -f $RD/*_part*_even.mrc $RD/*_part*_odd.mrc $RD/rho_*
    st $RP > $D/dump.txt; labels $D/dump.txt > $D/labels.txt
    $PY $T/p2_metrics.py $D/dump.txt --gt $S/gt_particles.txt --true $S/all.txt --ref0 $W/lab_0.txt --prev $W/lab_$((it - 1)).txt --nsmall 3 --tag "$ARM after it $it" 2>&1 | grep -v Warning
    cp -f $D/labels.txt $W/lab_$it.txt; cp -f $RP $W/proj_$((it + 1)).simple
    awk -v d=$RD 'NR == FNR {if ($1 > 0) c[$1]++; next} END {for (s = 1; s <= k; s++) printf "%d %s/recvol_state%02d_even_unfil.mrc %s/recvol_state%02d_odd_unfil.mrc %d\n", s, d, s, d, s, c[s]}' k=$(wc -l < $W/list_$it.txt) $D/labels.txt /dev/null > $W/list_$((it + 1)).txt
    awk -v d=$RD '{printf " vol%d=%s/recvol_state%02d.mrc", $1, d, $1}' $W/list_$((it + 1)).txt > $W/rawrefs_$((it + 1)).txt
  done
  F=$W/final; mkdir -p $F; cp -f $W/list_$((NIT + 1)).txt $F/units_list.txt
  [ -s $F/tree/merge_tree.txt ] || $PY $T/state_merge.py $F/tree --list $F/units_list.txt --link complete --ncut 2 > $F/tree.log 2>&1
  $PY $T/tree_wiener5.py $F/tw --list $F/units_list.txt --labels $W/lab_$NIT.txt --tree $F/tree/merge_tree.txt --gt $S/gt_particles.txt --truth $S/maps --mskdiam $MD 2>&1 | grep -v Warning | grep -v "smallest leaf" | sed "s/^/  $ARM final maps: /"
done
NR=$N/refined_states.txt; [ -s $NR ] && { RDN=$(ls -td $N/*_refine3D | head -1); st $RDN/r.simple > $O/dump_native.txt; $PY $T/p2_metrics.py $O/dump_native.txt --gt $S/gt_particles.txt --true $S/all.txt --ref0 $O/raw/lab_0.txt --nsmall 3 --tag "native node refinement (one refine3D call)" 2>&1 | grep -v Warning; }
log DONE
