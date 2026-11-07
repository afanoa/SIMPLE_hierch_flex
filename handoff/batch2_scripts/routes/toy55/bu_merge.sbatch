#!/bin/bash
#SBATCH -J klem_bumerge
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# generic merge-up along the bottom-up tree (bu_merge.py): levels G = N/2, N/4, ..., 1 of the leaf tree; per group the units are tested by
# splitbench8 on the even/odd halves of ONE stack (cand files carry labels on that half's rows only), no re-assignment, no blend pruning.
# env TAG (sim<TAG>/all.mrcs, lat_<TAG>/ptcl3D_r0.txt), LEAVES (per-row leaf labels), TREE (state_merge.py merge_tree.txt of the leaves),
# PROJ (project with the rows' eo), OUT, CTFABS (0), FOCUS (mask or empty), PAR (4), STRATA (1: null stratified by view x defocus, splitbench9), XEM (3: cross-half EM re-assignment drains states made of shared noise), ARCH (0); scorecard: GT, CLASSES, CODES (environment)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; TAG=${TAG:?}; O=${OUT:?}; mkdir -p $O; cd $O
log(){ echo "[$(date '+%F %T') bu-merge $TAG] $*"; }; S=${SIM:-$W/sim$TAG}; P=${POSES:-$W/lat_$TAG/ptcl3D_r0.txt}; CP=$W/benchcache_full_$TAG
[ -s $O/eo.txt ] || simple_exec prg=print_project_field projfile=${PROJ:?} oritype=ptcl3D mkdir=no 2>/dev/null | awk '{e = -1; for (i = 1; i <= NF; i++) if ($i ~ /^eo=/) {split($i, a, "="); e = a[2] + 0}; print e}' > $O/eo.txt
eo(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{e = -1; for (i = 1; i <= NF; i++) if ($i ~ /^eo=/) {split($i, a, "="); e = a[2] + 0}; print e}' > $O/eo.txt; }
if [ "$(sort -u $O/eo.txt | grep -c '^[01]$')" -lt 2 ]; then PR2=$(ls -td $(dirname $LEAVES)/*_refine3D_states 2>/dev/null | head -1)/p.simple; [ -s $PR2 ] && eo $PR2; fi   # input copies carry no split: the refinement's output project does
[ "$(sort -u $O/eo.txt | grep -c '^[01]$')" -eq 2 ] || { log "no even/odd split in the project rows: give PROJ = the refinement's output project"; exit 1; }
log "eo: $(awk '$1 == 0' $O/eo.txt | wc -l) even, $(awk '$1 == 1' $O/eo.txt | wc -l) odd of $(grep -c . $O/eo.txt) rows"
score(){ [ -n "${GT:-}" ] && $PY $W/scorecard.py $1 $GT ${CLASSES:+--classes $CLASSES} ${CODES:+--codes $CODES} --tag "$2" 2>&1 | grep -v Warning; }
one(){ d=$1; [ -s $d/children_B.txt ] && return 0
  $PY -u $W/splitbench9.py --sim $S,$S --poses $P,$P --cand $d/cand_A.txt,$d/cand_B.txt --out $d --cache $CP --ctfabs ${CTFABS:-0} --pair identity --xem ${XEM:-3} --softem 0 --arch ${ARCH:-0} --strata ${STRATA:-1} \
    $( [ -n "${FOCUS:-}" ] && echo "--focus $FOCUS" ) > $d/bench.log 2>&1 && echo "   $(basename $d): $(($(grep -c . $d/labels.txt) - 1)) units -> K* $(grep -a 'K\* =' $d/bench.log | tail -1 | grep -o '[0-9]*' | head -1)" \
    || echo "   $(basename $d) FAILED: $(grep -a -m1 'Error' $d/bench.log)"; }
U=$LEAVES; UL=-; N=$(awk '$1 > 0' $LEAVES | sort -n | uniq | wc -l); score $U "leaves"; g=$N; lv=0
while [ $g -gt 1 ]; do g=$(( (g + 1) / 2 )); lv=$((lv + 1)); L=$O/level$lv
  $PY $W/bu_merge.py prep $U $UL $O/eo.txt ${TREE:?} $L $g | sed "s/^/level $lv: /"
  first=1; while read d us; do case "$d" in \#*) continue;; esac; [ -s $d/cand_A.txt ] || continue
    if [ $first = 1 ]; then one $d; first=0; [ -s ${CP}_A_c48.npy ] && [ ! -L ${CP}_B_c48.npy ] && ln -sf ${CP}_A_c48.npy ${CP}_B_c48.npy.tmp && mv -f ${CP}_B_c48.npy.tmp ${CP}_B_c48.npy; continue; fi   # same stack: one cache
    while [ $(jobs -rp | wc -l) -ge ${PAR:-4} ]; do sleep 5; done; one $d &
  done < $L/nodes.txt; wait
  $PY $W/bu_merge.py collect $U $UL $L $L/out | sed "s/^/level $lv: /"; U=$L/out_units.txt; UL=$L/out_leaves.txt; score $U "level $lv (tree cut $g)"
done
cp -f $U $O/final_units.txt; cp -f $UL $O/final_leaves.txt; log "DONE: $(awk '$1 > 0' $U | sort -n | uniq | wc -l) states"
