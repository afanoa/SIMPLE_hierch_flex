#!/bin/bash
#SBATCH -p norm
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# Option (a) of the quality screen on labelled node_motion candidates, each node in its own frame (recB0/b.simple, the bins' box).
# env TREE (run dir), KIND (cdisc | igm | rnasep | other), NODES (comma list), NAME, EXCL (optional node:cand duplicates). Per node one
# quality_ref.sbatch reference (2 at a time), then qa_ratio.py cases, quality_screen3.py (and fa_maps10.py --dock for rnasep), qa_ratio.py
# table -> rnasep_klem/scorecard/quality_a/NAME/ratios.txt
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; FA=$K/full692_analysis
H=${TREE:?}; H=${H%/}; [ -d $H/1_hetero_tree ] && H=$H/1_hetero_tree; OUT=$K/scorecard/quality_a/${NAME:?}; RP=$K/scorecard/quality_ref/qa_${NAME}_
[ -e $OUT ] && mv $OUT ${OUT}_old_$(date +%s); mkdir -p $OUT; cd $OUT
log(){ echo "[$(date "+%F %T") qa_validate $NAME] $*"; }
for X in ${NODES//,/ }; do N=$H/node_$X; C=$(ls -d $N/motion/c[0-9]* | head -1)
  ( PROJ=$N/recB0/b.simple TAG=qa_${NAME}_$X RECSH=$C/hetero_step.sh MASK=$N/mask.mrc NP=2 NT=8 bash $T/quality_ref.sbatch > ref_$X.log 2>&1 ) &
  while [ $(jobs -rp | wc -l) -ge 2 ]; do sleep 5; done
done; wait
for X in ${NODES//,/ }; do log "reference $X: $(grep -a QS3 ${RP}$X/dBref.txt 2>/dev/null || echo FAILED, see ref_$X.log)"; done
$PY $T/qa_ratio.py cases $H ${KIND:?} $NODES $OUT $RP ${EXCL:-} > cases.log 2>&1 || { log "cases FAILED: $(tail -2 cases.log)"; exit 1; }
log "$(wc -l < cases.txt) candidates"; split -n l/4 -d cases.txt qcase_
for f in qcase_0[0-9]; do ( OMP_NUM_THREADS=4 $PY $T/quality_screen3.py $f < qscreen.txt > $f.qs3 2>&1 ) & done
if [ $KIND = rnasep ]; then split -n l/4 -d dock_maplist.txt dock_part_
  for f in dock_part_0[0-9]; do [ -s $f ] && ( OMP_NUM_THREADS=4 $PY -u $FA/fa_maps10.py $f dock_${f#dock_part_} --dock 2>&1 | grep -v "Warning\|warn(" > $f.log ) & done; fi
wait; cat qcase_0[0-9].qs3 > qs3.txt
if [ $KIND = rnasep ]; then mkdir -p dock; { head -1 dock_00/maps_table.txt; for d in dock_0[0-9]; do tail -n +2 $d/maps_table.txt 2>/dev/null; done; } > dock/maps_table.txt; fi
$PY $T/qa_ratio.py table $H $KIND $NODES $OUT $RP > table.log 2>&1 || log "table FAILED: $(tail -2 table.log)"; cat ratios.txt; log DONE
