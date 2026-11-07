#!/bin/bash
# masked split test of one f692e board (region board or zoom): usage prune_one.sh BOARD_DIR NAME OUTDIR
# (env FLOOR 500 particles per half; poses = the f692e family-refined poses mapped to the halves)
B=$1; NAME=$2; O=$3; mkdir -p $O
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; P=$K/dtree_rnasep692_f692e/poses_famrefined
CONS=/mnt/beegfs/elmlund/afan/rnasep_states/25_refine3D_states/rec_final_state01.mrc; CL=$W/maps160/docked_all/conf18_emd74209_m160.mrc
MASK=$(grep -h "^focus_mask" $B/board/board.txt 2>/dev/null | awk '{print $2}'); [ -n "$MASK" ] && [ -s "$MASK" ] || MASK=none
$PY $W/prune_prep.py $B/level1/states_clusters.txt $W/simrnasep692_A/meta_rows.txt $W/simrnasep692_B/meta_rows.txt ${FLOOR:-500} $O > $O/prep.log 2>&1 || { echo "PRUNE $NAME prep FAILED: $(tail -1 $O/prep.log)"; exit 0; }
nt=$(grep -vc "^#" $O/labels.txt)
if [ "$nt" -ge 2 ]; then
  $PY -u $W/splitbench8.py --sim $W/simrnasep692_A,$W/simrnasep692_B --poses ${P}_A.txt,${P}_B.txt --cand $O/cand_A.txt,$O/cand_B.txt --out $O \
    --cache $W/benchcache_rnasep692 --ctfabs 1 --pair identity --xem 0 --softem 0 $( [ $MASK != none ] && echo "--focus $MASK" ) > $O/bench.log 2>&1 || echo "PRUNE $NAME bench FAILED: $(grep -a -m2 'Error\|Traceback' $O/bench.log)"
else cp $O/cand_A.txt $O/children_A.txt; cp $O/cand_B.txt $O/children_B.txt; fi
$PY $W/arm_readout.py $CONS $CL $O/arm.txt $B/level1/states_dir > $O/arm.log 2>&1
$PY $W/prune_post.py $O $B/level1/states_dir $MASK $O/arm.txt $NAME 2>&1 | grep -v Warning
