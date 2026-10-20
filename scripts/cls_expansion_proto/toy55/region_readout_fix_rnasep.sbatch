#!/bin/bash
#SBATCH -J rnasep_region_readout
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Fix-up for the overnight region boards (the map copies had used an unpadded seq): waits for both drivers to finish,
# then for every region dir copies the level-1 and level-2 reconstructions under their proper names and redoes the
# identification (deposited_match.py) and the contact table (vis.py) per state; then the landscape count.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KF=$A/rnasep_klem/expand1/fam1_pipeline/klem_fb; D=$A/rnasep_klem/knobfree; VIS=$KF/closed_focal/lever/vis.py; log(){ echo "[$(date '+%F %T') rf] $*"; }
n=0; while [ $n -lt 300 ] && { ! grep -q "kf\] DONE" $A/rnasep_klem/rnasep_knobfree.53904039.out 2>/dev/null || ! grep -q "kf\] DONE" $A/rnasep_klem/rnasep_knobfree2.53904042.out 2>/dev/null; }; do sleep 60; n=$((n+1)); done
for R in $(ls -d $D/region_* | sort); do
  reg=$(basename $R); echo "===================== $reg"
  for L in $R/level1 $R/level2_s1 $R/level2_s2 $R/level2_s3; do
    [ -d $L/rec ] || continue; RR=$(ls -td $L/rec/*_reconstruct3D 2>/dev/null | head -1); [ -n "$RR" ] || { echo "  $(basename $L): no reconstruction"; continue; }
    mkdir -p $L/maps $L/states_dir; rm -f $L/states_dir/*.mrc; cnt=0
    for f in $RR/recvol_state??.mrc; do b=$(basename $f .mrc); [[ $b == *_even || $b == *_odd ]] && continue; cp -f $f $L/maps/; cp -f ${f%.mrc}_even.mrc ${f%.mrc}_odd.mrc $L/maps/ 2>/dev/null; cp -f $f $L/states_dir/state_${b#recvol_state}.mrc; cnt=$((cnt+1)); done
    echo "== $reg $(basename $L): $cnt maps; identification (deposited maps used only here)"; $PY $W/deposited_match.py $L/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc --save $L/cell_matches.txt 2>&1 | grep -v Warning | tail -13
    args=(); for f in $(ls $L/maps/recvol_state??.mrc 2>/dev/null); do s=$(basename $f .mrc); s=${s#recvol_state}; nsel=$(awk -v s=$((10#$s)) '$1==s' $L/states_clusters.txt 2>/dev/null | wc -l); [ "$(basename $L)" = level1 ] || nsel=$(awk '$1==1' $L/states_top.txt 2>/dev/null | wc -l); args+=("$f" "$reg $(basename $L) state $((10#$s)) ($nsel)"); done
    [ ${#args[@]} -gt 0 ] && $PY $VIS "${args[@]}" 2>&1 | grep -v Warning | grep -v "^  deposited"
  done
done
log "landscape count over all region boards (level 1)"; $PY $W/count_states.py $(ls -d $D/region_*/board) --sim 0.6,0.7,0.8 2>&1 | grep -v Warning | tail -8
log DONE
