#!/bin/bash
#SBATCH -J score_tree
#SBATCH -p norm
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# score_tree: the one-page usefulness scorecard of a hetero_tree run (score_tree.sh submits it). env TREE (the run dir holding
# 1_hetero_tree/), KIND (rnasep | cdisc | igm | cahra | other), PAIR (optional: a second seed's run dir). Reads the run only; writes
# rnasep_klem/scorecard/<run>/: jobs.txt (E, sacct), dump.npz (the root stage-B corr / dfx / poses), splits.txt (B, split_effect.py),
# leaves.txt (leaf_stats.py), fa_leaves/ + fa_split*/ + recall.txt + occupancy.txt (A, rnasep: fa_maps10.py --dock, fa_recall.py,
# fa_occ.py), truth.txt (C, toys), pair.txt (D), then scorecard.txt (scorecard_page.py). The blur flag's reference: quality_ref.sbatch on
# the root at the stage-B box (scorecard/quality_ref/<run>_root_b<box>/); a node-rows-only root dump is mapped onto project rows
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; FA=$K/full692_analysis
R=${TREE:?}; R=${R%/}; H=$R/1_hetero_tree; [ -d $H ] || { H=$R; R=$(dirname $H); }; NAME=$(basename $R); OUT=$K/scorecard/$NAME; mkdir -p $OUT; cd $OUT
log(){ echo "[$(date "+%F %T") score_tree $NAME] $*"; }
sacct -u $USER -S 2026-09-01 -X -n -P --format=JobID,JobName,WorkDir,ElapsedRaw,AllocCPUS,Start,End 2>/dev/null | awk -F'|' -v p="$R" 'index($3, p) == 1' > jobs.txt
log "$(wc -l < jobs.txt) jobs ran in $R"
( simple_exec prg=print_project_field projfile=$H/node_r/stageB/b.simple oritype=ptcl3D 2>/dev/null | $PY $T/dump_ptcl.py dump.npz > dump.log 2>&1 ) &
DUMP=$!
if [ -s $H/final_maps/units_index.txt ]; then bash $FA/fa_maplist_final.sh $H raw | awk '$2 != "node"' > leaves_maplist.txt; else : > leaves_maplist.txt; fi
$PY $T/split_maps.py $H > split_maplist.txt
if [ "${KIND:?}" = rnasep ]; then
  [ -s leaves_maplist.txt ] && ( OMP_NUM_THREADS=4 $PY -u $FA/fa_maps10.py leaves_maplist.txt fa_leaves --dock --title "$NAME leaves" 2>&1 | grep -v "Warning\|warn(" > fa_leaves.log
    $PY $FA/fa_recall.py fa_leaves > recall.txt 2>&1; $PY $FA/fa_occ.py fa_leaves > occupancy.txt 2>&1 ) &
  split -n l/8 -d split_maplist.txt split_part_; SPIDS=()
  for f in split_part_0*; do [ -s $f ] || continue; ( OMP_NUM_THREADS=2 $PY -u $FA/fa_maps10.py $f fa_${f#split_part_} --dock 2>&1 | grep -v "Warning\|warn(" > fa_${f#split_part_}.log ) & SPIDS+=($!); done
fi
case $KIND in
  cdisc) ( $PY $T/validate_tree3.py $NAME $H $K/factorial/junk_test/simsites_truth_jt.txt sites > truth.txt 2>&1 ) & ;;
  igm)   ( $PY $T/validate_tree3.py $NAME $H $K/factorial/armsigm1_truth/arms_truth.txt arms > truth.txt 2>&1
           $PY $T/arm_score.py $H $K/factorial/armsigm1_truth/arms_truth.txt > armscore.txt 2>&1 ) & ;;
  cahra) ( $PY $T/validate_conf.py $NAME $H /mnt/beegfs/elmlund/EMPIAR_downloads/13991/data/source.star > truth.txt 2>&1 ) & ;;
  qtoy)  ( $PY $T/qtoy_score.py $H ${TRUTH:-$T/simqtoy/truth_qtoy.txt} > qtoyscore.txt 2>&1 ) & ;;
esac
[ -n "${PAIR:-}" ] && ( $PY /mnt/beegfs/elmlund/afan/ole_12707/scripts/tree_pair_verify.py $R $PAIR > pair.txt 2>&1 ) &
wait $DUMP
$PY - $H <<'EOF' >> dump.log 2>&1
import sys, numpy as np
H = sys.argv[1]; D = dict(np.load('dump.npz')); r0 = np.loadtxt(H + '/node_r/rows.txt', dtype=int); nP = sum(1 for _ in open(H + '/hetero_tree_labels.txt')); nd = len(D['state'])
if nd != nP:
    assert nd == len(r0) and r0.max() <= nP, 'dump %d rows: neither the %d project rows nor the %d root rows' % (nd, nP, len(r0))
    for k in D: v = np.full(nP, np.nan); v[r0 - 1] = D[k]; D[k] = v
    np.savez('dump.npz', **D); print('%d root rows mapped onto %d project rows through node_r/rows.txt' % (nd, nP))
open('dump_mode.txt', 'w').write('rows\n' if nd != nP else 'full\n')
EOF
log "dump: $(tail -1 dump.log)"
B88=$(od -An -t d4 -N4 $H/node_r/stageB/recvol_state01.mrc | tr -d ' '); MB=$(od -An -t d4 -N4 $H/node_r/mask.mrc | tr -d ' '); QT=${NAME}_root_b$B88
( PROJ=$H/node_r/stageB/b.simple TAG=$QT RECSH=$H/node_r/recB0/hetero_step.sh MASK=$H/node_r/mask.mrc CROP=$([ "$B88" != "$MB" ] && echo $B88) NP=4 NT=8 bash $T/quality_ref.sbatch > qref.log 2>&1 ) &
QREF=$!
[ -s leaves_maplist.txt ] && ( $PY $T/leaf_stats.py $H dump.npz leaves_maplist.txt leaves.txt > leaf_stats.log 2>&1 ) &
if [ $KIND = rnasep ]; then [ ${#SPIDS[@]} -gt 0 ] && wait ${SPIDS[@]}
  { head -1 fa_00/maps_table.txt; for d in fa_0*; do tail -n +2 $d/maps_table.txt 2>/dev/null; done; } > split_maps_table.txt; FAOPT="--fa split_maps_table.txt"; else FAOPT=""; fi
wait $QREF; QR=$K/scorecard/quality_ref/$QT/dBref.txt; log "root quality reference: $(grep -a QS3 $QR 2>/dev/null || echo FAILED, see qref.log)"
$PY $T/split_effect.py $H dump.npz splits.txt $FAOPT --dbref $QR > split_effect.log 2>&1
wait
$PY $T/scorecard_page.py $OUT $H $NAME $KIND > /dev/null 2>&1; log "DONE: $OUT/scorecard.txt ($(wc -l < scorecard.txt) lines)"
