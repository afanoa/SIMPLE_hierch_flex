#!/bin/bash
# one zoom: a board on the particles of level-1 state S of region RB (family FAM of run R), its supported sub-states,
# each reconstructed (local parts). args: R FAM RB S ; env: W PY NTHR (12) NPZ (2 rec parts) SIM (0.7) ZMAX (10000)
set -uo pipefail
R=$1; FAM=$2; RB=$3; S=$4; W=${W:-/mnt/beegfs/elmlund/afan/routes/toy55}; PY=${PY:-/mnt/nasapps/production/python/3.12/bin/python3}
D=$R/fam$FAM; REG=$(basename $RB | sed 's/^region//'); Z=$D/zoomregion${REG}_s$(printf %02d $S); L=$Z/level1; mkdir -p $L; log(){ echo "[$(date '+%F %T') zoom f$FAM/$REG/s$S] $*"; }
[ -f $Z/DONE ] && exit 0
PS=$(ls -td $D/refine/*_refine3D_states 2>/dev/null | head -1)/p.simple; [ -f "$PS" ] || PS=${SRC:?}; META=$D/meta_fb.npz; NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))")
MK=$(ls $D/regions_motion/region_${REG#motion_}.mrc 2>/dev/null); case $REG in motion_*) MK=$D/regions_motion/region_${REG#motion_}.mrc;; *) MK=$D/regions/region${REG}.mrc;; esac
awk -v s=$S '$1==s{print NR}' $RB/level1/states_clusters.txt > $Z/subset.txt; n=$(grep -c . $Z/subset.txt); [ $n -ge 2000 ] || { log "$n particles, no zoom"; touch $Z/DONE; exit 0; }
log "zoom board on $n particles (focus $(basename $MK))"
[ -s $Z/board/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=${NTHR:-12} N1=7 N2=7 NFREE=1 CTFABS=${CTFABS:-1} PSCALE=${PSCALE:-1} FOCUSMASK=$MK META=$META bash $W/fboard.sh $Z/board $PS $Z/subset.txt "${ZLPS:-12 10 10 8 8 8}" "${ZSIG:-1.2 1.0 1.0 0.8 0.7 0.7}" > $Z/board.out 2>&1 || { log "board FAILED: $(tail -2 $Z/board.out | cut -c1-120)"; touch $Z/DONE; exit 1; }
$PY $W/node_states2.py $Z/board $L --nrows $NROWS --sim ${SIM:-0.7} --nmax ${ZMAX:-10000} --support --maxk 10 > $Z/states.out 2>&1; K=$(grep -vc "^#" $L/clusters.txt 2>/dev/null || echo 0)
[ "$K" -ge 1 ] || { log "no supported sub-states"; touch $Z/DONE; exit 0; }
log "$K sub-states: $(grep -v '^#' $L/clusters.txt | awk '{printf "%s(%s) ", $4, $1}')"
mkdir -p $L/rec $L/maps $L/states_dir && cd $L/rec && cp -f $PS p.simple && simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=${NPZ:-2} nthr=${NTHR:-12} mkdir=no > up.log 2>&1
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$L/states_clusters.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; touch $Z/DONE; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$K nparts=${NPZ:-2} nthr=${NTHR:-12} > rec.log 2>&1 || { log "reconstruct3D FAILED"; touch $Z/DONE; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); for s in $(seq -f "%02g" 1 $K); do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $L/maps/ 2>/dev/null; done; cp -f $RR/recvol_state$s.mrc $L/states_dir/state_$s.mrc; done
rm -rf $RR/*_part*.mrc 2>/dev/null; touch $Z/DONE; log "DONE ($K sub-states)"
