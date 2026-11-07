#!/bin/bash
# f692 regions stage (truth-free node recipe with the Fortran board) for one family: regions of the feedback board
# (spatial components + motion-pattern groups), one focal board per region at the refined poses, states = clusters of
# reproducible cells (node_states.py), top particles of every state reconstructed at the refined poses, level 2 for
# the largest states (within-state board with the inherited region, ranked along the node's own main axis). The
# deposited maps enter only in the identification readouts.
# env: R, FAM, META (meta_fb of the family), MAXREG (6), NL2 (3), NMAX (30000), N2 (10000), SIM (0.7), NPREC (4)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
D=$R/fam$FAM; SUB=$R/expand/family${FAM}_subset.txt; [ -s $R/maxreg_fam$FAM.txt ] && MAXREG=$(cat $R/maxreg_fam$FAM.txt); PS=$(ls -td $D/refine/*_refine3D_states | head -1)/p.simple; log(){ echo "[$(date '+%F %T') reg$FAM] $*"; }
NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))"); LPSR="12 12 10 10 8 8 8 8 8"; SIGR="1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7"
rec(){ dir=$1; sf=$2; ns=$3; mkdir -p $dir/rec && cd $dir/rec && cp -f $PS p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$sf mkdir=no > sel.log 2>&1 || { log "selection FAILED ($dir)"; cd $D; return 1; }
  simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=${NPREC:-4} nthr=8 > rec.log 2>&1 || { log "reconstruct3D failed once ($dir), retrying"; sleep 30; simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=${NPREC:-4} nthr=8 > rec2.log 2>&1 || { log "reconstruct3D FAILED twice ($dir): $(grep -a -m2 ERROR rec2.log | cut -c1-200)"; cd $D; return 1; }; }
  RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $dir/maps $dir/states_dir; for s in $(seq -f "%02g" 1 $ns); do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $dir/maps/ 2>/dev/null; done; cp -f $RR/recvol_state$s.mrc $dir/states_dir/state_$s.mrc; done; cd $D; }
# ---- regions of the feedback board
[ -s $D/regions/region_1.mrc ] || $PY $W/derive_regions.py $D/klem_fb $D/regions --maxreg 2 2>&1 | grep -v Warning | tail -4
[ -s $D/regions_motion/region_1.mrc ] || $PY $W/derive_regions2.py $D/klem_fb $D/regions_motion --k 5 2>&1 | grep -v Warning | tail -7
MASKS=$( { ls $D/regions/region_?.mrc 2>/dev/null; ls $D/regions_motion/region_?.mrc 2>/dev/null; } | head -${MAXREG:-6}); log "$(echo $MASKS | wc -w) region masks: $(for m in $MASKS; do echo -n "$(basename $(dirname $m))/$(basename $m .mrc) "; done)"
export OMP_NUM_THREADS=8
for MK in $MASKS; do
  REG=$(basename $(dirname $MK) | sed 's/regions//; s/^_//')$(basename $MK .mrc | sed 's/region//'); RB=$D/region$REG; mkdir -p $RB
  log "(1) region $REG: focal board at the refined poses"; [ -s $RB/board/klem_latent.txt ] || FOCUSMASK=$MK META=$META bash $W/fboard.sh $RB/board $PS $SUB "$LPSR" "$SIGR" || { log "board FAILED ($REG)"; continue; }
  log "(2) states from the reproducible cell clusters"; $PY $W/node_states.py $RB/board $RB/level1 --nrows $NROWS --sim ${SIM:-0.7} --nmax ${NMAX:-30000} 2>&1 | grep -v Warning | tail -3
  KS=$(grep -vc "^#" $RB/level1/clusters.txt); [ "$KS" -ge 1 ] || { log "no states ($REG)"; continue; }
  [ -s $RB/level1/maps/recvol_state01.mrc ] || rec $RB/level1 $RB/level1/states_clusters.txt $KS || continue
  echo "== identification of the level-1 states of region $REG (readout)"; $PY $W/deposited_match.py $RB/level1/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc --save $RB/level1/state_matches.txt 2>&1 | grep -v Warning | tail -$((KS + 2))
  log "(3) level 2 for the largest ${NL2:-3} states (inherited region, own main axis)"
  for s in $(seq 1 ${NL2:-3}); do [ $s -le $KS ] || break; L=$RB/level2_s$s; mkdir -p $L; awk -v s=$s '$1==s{print NR}' $RB/level1/states_clusters.txt > $L/subset.txt; n=$(grep -c . $L/subset.txt); [ $n -ge 4000 ] || { log "state $s: $n particles, no level 2"; continue; }
    [ -s $L/board/klem_latent.txt ] || FOCUSMASK=$MK META=$META bash $W/fboard.sh $L/board $PS $L/subset.txt "$LPSR" "$SIGR" || continue
    n2=$(( n * 4 / 10 )); [ $n2 -le ${N2:-10000} ] || n2=${N2:-10000}; $PY $W/rank_by_target.py $L/board --target pc1 --n $n2 --nrows $NROWS --out $L 2>&1 | grep -v Warning | tail -2
    [ -s $L/maps/recvol_state01.mrc ] || rec $L $L/states_top.txt 1 || continue
    echo "== identification of region $REG state $s level 2 (top $n2, readout)"; $PY $W/deposited_match.py $L/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -3
  done
done
log DONE
