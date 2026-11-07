#!/bin/bash
#SBATCH -J rnasep_f692_region
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# ONE region of the truth-free node recipe, run as its own job so the regions of a family run in parallel: focal board
# at the refined poses with the region as focus, states = clusters of reproducible cells (node_states.py), the top
# particles of every state reconstructed at the refined poses, level 2 for the largest states (inherited region, ranked
# along the node's own main axis). Writes $RB/DONE at the end. env: R, FAM, MK (mask), REG (tag), META (meta_fb),
# FB (build), NL2 (3), NMAX (30000), N2 (10000), SIM (0.7), NPREC (4), NTHR (32)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
D=$R/fam$FAM; SUB=$R/expand/family${FAM}_subset.txt; PS=$(ls -td $D/refine/*_refine3D_states | head -1)/p.simple; RB=$D/region$REG; mkdir -p $RB; log(){ echo "[$(date '+%F %T') reg$FAM/$REG] $*"; }
NTOP2=${NTOP2:-${N2:-10000}}; unset N2   # the level-2 particle count; N2 must never reach fboard.sh (its lattice dimension)
NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))"); LPSR="12 12 10 10 8 8 8 8 8"; SIGR="1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7"; export NTHR=${NTHR:-16}
rec(){ dir=$1; sf=$2; ns=$3; mkdir -p $dir/rec && cd $dir/rec && cp -f $PS p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$sf mkdir=no > sel.log 2>&1 || { log "selection FAILED ($dir)"; cd $RB; return 1; }
  simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=${NPREC:-4} nthr=8 > rec.log 2>&1 || { log "reconstruct3D failed once ($dir), retrying"; sleep 30; simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=${NPREC:-4} nthr=8 > rec2.log 2>&1 || { log "reconstruct3D FAILED twice ($dir)"; cd $RB; return 1; }; }
  RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $dir/maps $dir/states_dir; for s in $(seq -f "%02g" 1 $ns); do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $dir/maps/ 2>/dev/null; done; cp -f $RR/recvol_state$s.mrc $dir/states_dir/state_$s.mrc; done; cd $RB; }
export OMP_NUM_THREADS=8
log "(1) focal board at the refined poses, mask $MK"; [ -s $RB/board/klem_latent.txt ] || FOCUSMASK=$MK META=$META bash $W/fboard.sh $RB/board $PS $SUB "$LPSR" "$SIGR" || { log "board FAILED"; exit 1; }
log "(2) states from the reproducible cell clusters"; $PY $W/node_states2.py $RB/board $RB/level1 --nrows $NROWS --sim ${SIM:-0.7} --nmax ${NMAX:-30000} --maxk ${MAXK:-30} --support 2>&1 | grep -v Warning | tail -3
KS=$(grep -vc "^#" $RB/level1/clusters.txt); [ "$KS" -ge 1 ] || { log "no states"; touch $RB/DONE; exit 0; }
[ "${NOREC:-0}" = 1 ] && { log "NOREC: $KS states from the board, no reconstructions or level 2"; touch $RB/DONE; exit 0; }
if [ "${NOREC_L1:-0}" = 1 ]; then log "NOREC_L1: level-1 states kept as lists only (no reconstructions, no identification), level 2 continues"; else
[ -s $RB/level1/maps/recvol_state01.mrc ] || rec $RB/level1 $RB/level1/states_clusters.txt $KS || { touch $RB/DONE; exit 1; }
echo "== identification of the level-1 states of region $REG (readout)"; $PY $W/deposited_match.py $RB/level1/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc --save $RB/level1/state_matches.txt 2>&1 | grep -v Warning | tail -$((KS + 2))
fi
log "(3) level 2 for the largest ${NL2:-3} states (inherited region, own main axis)"
for s in $(seq 1 ${NL2:-3}); do [ $s -le $KS ] || break; L=$RB/level2_s$s; mkdir -p $L; awk -v s=$s '$1==s{print NR}' $RB/level1/states_clusters.txt > $L/subset.txt; n=$(grep -c . $L/subset.txt); [ $n -ge 4000 ] || { log "state $s: $n particles, no level 2"; continue; }
  [ -s $L/board/klem_latent.txt ] || NPARTS= FOCUSMASK=$MK META=$META bash $W/fboard.sh $L/board $PS $L/subset.txt "$LPSR" "$SIGR" || continue   # level 2 on the host threads: the distributed board has a ~90 s/iteration exchange floor
  n2=$(( n * 4 / 10 )); [ $n2 -le $NTOP2 ] || n2=$NTOP2; $PY $W/rank_by_target.py $L/board --target pc1 --n $n2 --nrows $NROWS --out $L 2>&1 | grep -v Warning | tail -2
  [ -s $L/maps/recvol_state01.mrc ] || rec $L $L/states_top.txt 1 || continue
  echo "== identification of region $REG state $s level 2 (top $n2, readout)"; $PY $W/deposited_match.py $L/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -3
done
touch $RB/DONE; log DONE
