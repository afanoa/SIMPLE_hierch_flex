#!/bin/bash
# f692 regions stage, INLINE (LOCALDEV): the regions of a family run one after the other inside the family job, derive the family's regions from its feedback board, submit one f692_region_one job
# per region (up to MAXREG at a time, PARMAX concurrent), wait for their DONE markers. Same outputs and layout as the
# sequential f692_regions.sh. env: R, FAM, META, MAXREG (6), PARMAX (3), NL2, NMAX, N2, SIM, NPREC, FB, NTHR
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55
D=$R/fam$FAM; [ -s $R/maxreg_fam$FAM.txt ] && MAXREG=$(cat $R/maxreg_fam$FAM.txt); [ -s $R/parmax.txt ] && PARMAX=$(cat $R/parmax.txt); [ -s $R/nthr.txt ] && NTHR=$(cat $R/nthr.txt); log(){ echo "[$(date '+%F %T') regpar$FAM] $*"; }
[ -s $D/regions/region_1.mrc ] || $PY $W/derive_regions.py $D/klem_fb $D/regions --maxreg 2 2>&1 | grep -v Warning | tail -4
[ -s $D/regions_motion/region_1.mrc ] || $PY $W/derive_regions2.py $D/klem_fb $D/regions_motion --k 5 2>&1 | grep -v Warning | tail -7
if [ "${SOLVENT:-0}" = 1 ]; then [ -s $D/solvent_mask.mrc ] || $PY $W/derive_solvent.py $D/klem_fb/cells $D/solvent_mask.mrc 2>&1 | grep -v Warning | tail -1; export SOLVENTMASK=$D/solvent_mask.mrc; log "region boards with the family solvent mask $SOLVENTMASK"; else unset SOLVENTMASK; fi   # SOLVENT=1: solvent mask from the feedback cells on the region boards
MASKS=$( { ls $D/regions/region_?.mrc 2>/dev/null; ls $D/regions_motion/region_?.mrc 2>/dev/null; } | head -${MAXREG:-6}); log "$(echo $MASKS | wc -w) region masks (up to ${PARMAX:-3} jobs at a time)"
for MK in $MASKS; do
  REG=$(basename $(dirname $MK) | sed 's/regions//; s/^_//')$(basename $MK .mrc | sed 's/region//'); RB=$D/region$REG; mkdir -p $RB
  [ -f $RB/DONE ] && { log "region $REG already done"; continue; }
  log "region $REG inline"; R=$R FAM=$FAM MK=$MK REG=$REG META=$META FB=${FB:-$R/klem2_bin} NL2=${NL2:-3} NMAX=${NMAX:-30000} NTOP2=${N2:-10000} SIM=${SIM:-0.7} NPREC=${NPREC:-4} NTHR=${NTHR:-48} LOCALDEV=1 bash $W/f692_region_one3.sbatch 2>&1 | grep -v Warning
done
log "regions finished: $(for MK in $MASKS; do REG=$(basename $(dirname $MK) | sed 's/regions//; s/^_//')$(basename $MK .mrc | sed 's/region//'); [ -f $D/region$REG/DONE ] && echo -n "region$REG "; done)"
