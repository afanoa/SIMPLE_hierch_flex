#!/bin/bash
#SBATCH -J klem_cheapkids_real
#SBATCH -p norm
#SBATCH -c 64
#SBATCH --mem=240G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# cheap leaf route on REAL rflex nodes (cheapkids_rf.py: node flex latent centred per direction bin, view-balanced quantile bins along PC1,
# PC2 and a 4x4 grid of the two) -> reconstruct3D at the parent's poses -> closed-arm readout per child. Maps (full) in
# view_rflex/cheap_<RUN>_<NODE>_<axis>/ (lp 10). env RUN (rflex_f692rf), NODES ("r22 r21 r2"), K (10), G (4), MSKDIAM (232)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RUN=${RUN:-rflex_f692rf}; RD=$K_/$RUN
K=${K:-10}; G=${G:-4}; MD=${MSKDIAM:-232}; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); log(){ echo "[$(date '+%F %T') cheapkids-real $RUN] $*"; }
for NODE in ${NODES:-r22 r21 r2}; do
  O=$K_/cheapkids_${RUN#rflex_}_$NODE; mkdir -p $O && cd $O
  set -- $(awk -v u=$NODE '$1 == u' $RD/queue.txt); PP=${3:-}; ROWS=${4:-}; FLEX=$(ls -td $RD/node_$NODE/*_flex_pca | head -1)/flex_pca_coordinates.txt
  [ -s "$PP" ] && [ -s "$ROWS" ] && [ -s "$FLEX" ] || { log "node $NODE inputs missing"; continue; }
  [ -s $O/ptcl3D_parent.txt ] || simple_exec prg=print_project_field projfile=$PP oritype=ptcl3D mkdir=no > $O/ptcl3D_parent.txt 2> /dev/null
  for ax in 1 2 grid; do kk=$K; [ $ax = grid ] && kk=$G
    [ -s $O/children_$ax.txt ] || $PY $T/cheapkids_rf.py --oris $O/ptcl3D_parent.txt --rows $ROWS --flex $FLEX --out $O --k $kk --axis $ax 2>&1 | grep -v Warning | sed "s/^/node $NODE /"
    NS=$(awk '$1 > 0' $O/children_$ax.txt | sort -n | uniq | wc -l)
    if [ -z "$(ls $O/rec_$ax/*_reconstruct3D/recvol_state$(printf %02d $NS).mrc 2>/dev/null)" ]; then mkdir -p $O/rec_$ax && cd $O/rec_$ax
      cp -f $PP r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/children_$ax.txt mkdir=no > rsel.log 2>&1 \
        && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > rup.log 2>&1 \
        && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=$NS nparts=4 nthr=16 > rec.log 2>&1 || { log "node $NODE axis $ax reconstruct FAILED"; cd $O; continue; }
      RR=$(ls -td $O/rec_$ax/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*; cd $O; fi
    RR=$(ls -td $O/rec_$ax/*_reconstruct3D | head -1); V=$K_/view_rflex/cheap_${RUN#rflex_}_${NODE}_$ax
    : > $O/spec_$ax.txt; for s in $(seq 1 $NS); do b=$RR/recvol_state$(printf %02d $s); echo "child_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec_$ax.txt; done
    $PY $T/make_view.py $V $O/spec_$ax.txt > /dev/null 2>&1
    $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_$ax.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | tail -1 | sed "s/^/node $NODE axis $ax: /"
    awk -v ax=$ax 'NR == FNR {n[$1]++; next} !/^#/ {k = $1; sub(/.*child_/, "", k); sub(/\.mrc/, "", k); printf "   child %s (%d ptcls): closed-arm %.2f open-arm %.2f\n", k, n[k + 0], $2, $3}' $O/children_$ax.txt $O/arm_$ax.txt
  done
done
log DONE
