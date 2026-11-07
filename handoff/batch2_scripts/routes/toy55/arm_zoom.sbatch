#!/bin/bash
#SBATCH -J arm_zoom
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# arm-focused zoom on ONE half of the real gold-standard tree (env H = A|B, TAG, SEL = node:<tree node> | state:<final state id>):
# klem board on the selected particles focused on the arm-travel mask (deposited closed conformer 18 vs consensus,
# make_arm_mask.py), 7x7 + 1 free, 12 -> 8 A, |CTF|, pscale; supported sub-states (node_states2.py --support), each
# reconstructed with SIMPLE from this half's project; closed-arm readout (arm_readout.py) of every sub-state map.
# Run once per half: a closed sub-state counts only if both halves produce one.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local TOY_SMPD=1.464 OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; H=${H:?}; TAG=${TAG:-gs1}; SEL=${SEL:?}; MODE=${MODE:-rnasep692}
R=$K/synth_${MODE}_$H; PS=$R/src.simple; META=$R/meta_toy.npz; NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))")
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; CLOSED=$W/maps160/docked_all/conf18_emd74209_m160.mrc; T=$K/dtree_${MODE}_$TAG
Z=$T/armzoom_$(echo $SEL | tr ':.' '__')_$H; mkdir -p $Z; log(){ echo "[$(date '+%F %T') armzoom $SEL $H] $*"; }
[ -s $T/arm_mask.mrc ] || $PY $W/make_arm_mask.py $CONS $CLOSED $T/arm_mask.mrc 2>&1 | grep -v Warning
case $SEL in node:*) cp -f $R/$TAG/${SEL#node:}/subset.txt $Z/subset.txt;; state:*) awk -v s=${SEL#state:} '$1 == s {print NR}' $T/${FIN:-final_xem}/children_$H.txt > $Z/subset.txt;; esac
n=$(grep -c . $Z/subset.txt); log "$n particles, focus $T/arm_mask.mrc"; [ $n -ge 2000 ] || { log "too few particles"; exit 1; }
[ -s $Z/board/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=16 N1=${N1:-7} N2=${N2:-7} NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$T/arm_mask.mrc META=$META SEED=${SEED:-1} \
  bash $W/fboard.sh $Z/board $PS $Z/subset.txt "${ZLPS:-12 10 10 8 8 8}" "${ZSIG:-1.2 1.0 1.0 0.8 0.7 0.7}" > $Z/board.out 2>&1 || { log "board FAILED: $(tail -2 $Z/board.out | cut -c1-200)"; exit 1; }
$PY $W/node_states2.py $Z/board $Z/level1 --nrows $NROWS --sim ${SIM:-0.7} --nmax ${ZMAX:-10000} --support --maxk 10 > $Z/states.out 2>&1; KS=$(grep -vc "^#" $Z/level1/clusters.txt 2>/dev/null || echo 0)
log "$KS supported sub-states: $(grep -v '^#' $Z/level1/clusters.txt 2>/dev/null | awk '{printf "%s(%s) ", $4, $1}')"; [ "$KS" -ge 1 ] || exit 0
mkdir -p $Z/level1/rec $Z/level1/states_dir && cd $Z/level1/rec && cp -f $PS p.simple && simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=1 nthr=16 mkdir=no > up.log 2>&1
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$Z/level1/states_clusters.txt mkdir=no > sel.log 2>&1 && \
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$KS nparts=1 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 'ERROR\|STOP' rec.log sel.log | cut -c1-160)"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); for s in $(seq -f "%02g" 1 $KS); do cp -f $RR/recvol_state$s.mrc $Z/level1/states_dir/state_$s.mrc; done
$PY $W/arm_readout.py $CONS $CLOSED $Z/arm_readout.txt $Z/level1/states_dir 2>&1 | grep -v Warning
log DONE
