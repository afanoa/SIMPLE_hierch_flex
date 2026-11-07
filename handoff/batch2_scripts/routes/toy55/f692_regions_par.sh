#!/bin/bash
# f692 regions stage, PARALLEL: derive the family's regions from its feedback board, submit one f692_region_one job
# per region (up to MAXREG at a time, PARMAX concurrent), wait for their DONE markers. Same outputs and layout as the
# sequential f692_regions.sh. env: R, FAM, META, MAXREG (6), PARMAX (3), NL2, NMAX, N2, SIM, NPREC, FB, NTHR
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
D=$R/fam$FAM; [ -s $R/maxreg_fam$FAM.txt ] && MAXREG=$(cat $R/maxreg_fam$FAM.txt); [ -s $R/parmax.txt ] && PARMAX=$(cat $R/parmax.txt); [ -s $R/nthr.txt ] && NTHR=$(cat $R/nthr.txt); log(){ echo "[$(date '+%F %T') regpar$FAM] $*"; }
[ -s $D/regions/region_1.mrc ] || $PY $W/derive_regions.py $D/klem_fb $D/regions --maxreg 2 2>&1 | grep -v Warning | tail -4
[ -s $D/regions_motion/region_1.mrc ] || $PY $W/derive_regions2.py $D/klem_fb $D/regions_motion --k 5 2>&1 | grep -v Warning | tail -7
MASKS=$( { ls $D/regions/region_?.mrc 2>/dev/null; ls $D/regions_motion/region_?.mrc 2>/dev/null; } | head -${MAXREG:-6}); log "$(echo $MASKS | wc -w) region masks (up to ${PARMAX:-3} jobs at a time)"
JOBS=""; RBS=""
for MK in $MASKS; do
  REG=$(basename $(dirname $MK) | sed 's/regions//; s/^_//')$(basename $MK .mrc | sed 's/region//'); RB=$D/region$REG; mkdir -p $RB
  [ -f $RB/DONE ] && { log "region $REG already done"; continue; }
  while [ $(echo $JOBS | wc -w) -ge ${PARMAX:-3} ]; do   # concurrency cap: wait for one running region job to finish
    sleep 60; for j in $JOBS; do squeue -h -j $j > /dev/null 2>&1 && squeue -h -j $j | grep -q . || JOBS=$(echo $JOBS | tr ' ' '\n' | grep -vx $j | tr '\n' ' '); done
  done
  jid=$(sbatch --parsable -c ${NTHR:-32} --export=ALL,R=$R,FAM=$FAM,MK=$MK,REG=$REG,META=$META,FB=${FB:-$R/klem2_bin},NL2=${NL2:-3},NMAX=${NMAX:-30000},NTOP2=${N2:-10000},SIM=${SIM:-0.7},NPREC=${NPREC:-4},NTHR=${NTHR:-32},NPARTS=${NPARTS:-} -J rnasep_f692_reg${FAM}_$REG $W/f692_region_one.sbatch) || { log "sbatch FAILED for region $REG"; continue; }
  JOBS="$JOBS $jid"; RBS="$RBS $RB"; log "region $REG -> job $jid"
done
n=0; while [ $n -lt 720 ]; do left=0; for RB in $RBS; do [ -f $RB/DONE ] || left=$((left+1)); done; [ $left -eq 0 ] && break; sleep 60; n=$((n+1)); done
log "regions finished: $(for RB in $RBS; do [ -f $RB/DONE ] && echo -n "$(basename $RB) "; done)"
