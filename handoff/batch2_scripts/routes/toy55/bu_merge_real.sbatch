#!/bin/bash
#SBATCH -J klem_bumerge
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# bu_merge_real.sbatch = bu_merge.sbatch on REAL half-stacks (HALVES_A / HALVES_B with meta_rows.txt, per-half pose tables POSES_A/B,
# their image caches, |CTF|): candidates and children live in each half's row space and map back through the row maps.
# generic merge-up along the bottom-up tree (bu_merge.py): levels G = N/2, N/4, ..., 1 of the leaf tree; per group the units are tested by
# splitbench8 on the even/odd halves of ONE stack (cand files carry labels on that half's rows only), no re-assignment, no blend pruning.
# env TAG (sim<TAG>/all.mrcs, lat_<TAG>/ptcl3D_r0.txt), LEAVES (per-row leaf labels), TREE (state_merge.py merge_tree.txt of the leaves),
# PROJ (project with the rows' eo), OUT, CTFABS (0), FOCUS (mask or empty), PAR (4), STRATA (1: null stratified by view x defocus, splitbench9), XEM (3: cross-half EM re-assignment drains states made of shared noise), ARCH (0); scorecard: GT, CLASSES, CODES (environment)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; TAG=${TAG:?}; O=${OUT:?}; mkdir -p $O; cd $O
log(){ echo "[$(date '+%F %T') bu-merge-real $TAG] $*"; }; H=${HALVES:-$W/simrnasep692}; MRA=${H}_A/meta_rows.txt; MRB=${H}_B/meta_rows.txt; CP=${CACHE:-$W/benchcache_rnasep692}
log "half-stacks ${H}_A / ${H}_B: $(grep -c . $MRA) / $(grep -c . $MRB) rows"
score(){ [ -n "${GT:-}" ] && $PY $W/scorecard.py $1 $GT ${CLASSES:+--classes $CLASSES} ${CODES:+--codes $CODES} --tag "$2" 2>&1 | grep -v Warning; }
one(){ d=$1; [ -s $d/children_B.txt ] && return 0
  $PY -u $W/splitbench9.py --sim ${H}_A,${H}_B --poses ${POSES_A:?},${POSES_B:?} --cand $d/cand_A.txt,$d/cand_B.txt --out $d --cache $CP --ctfabs ${CTFABS:-1} --pair identity --xem ${XEM:-3} --softem 0 --arch ${ARCH:-0} --strata ${STRATA:-1} \
    $( [ -n "${FOCUS:-}" ] && echo "--focus $FOCUS" ) > $d/bench.log 2>&1 && echo "   $(basename $d): $(($(grep -c . $d/labels.txt) - 1)) units -> K* $(grep -a 'K\* =' $d/bench.log | tail -1 | grep -o '[0-9]*' | head -1)" \
    || echo "   $(basename $d) FAILED: $(grep -a -m1 'Error' $d/bench.log)"; }
U=$LEAVES; UL=-; N=$(awk '$1 > 0' $LEAVES | sort -n | uniq | wc -l); score $U "leaves"; g=$N; lv=0
while [ $g -gt 1 ]; do g=$(( (g + 1) / 2 )); lv=$((lv + 1)); L=$O/level$lv
  $PY $W/bu_merge.py prep $U $UL - ${TREE:?} $L $g $MRA $MRB | sed "s/^/level $lv: /"
  first=1; while read d us; do case "$d" in \#*) continue;; esac; [ -s $d/cand_A.txt ] || continue
    if [ $first = 1 ]; then one $d; first=0; continue; fi   # the first test alone (caches exist per half)
    while [ $(jobs -rp | wc -l) -ge ${PAR:-4} ]; do sleep 5; done; one $d &
  done < $L/nodes.txt; wait
  $PY $W/bu_merge.py collect $U $UL $L $L/out $MRA $MRB | sed "s/^/level $lv: /"; U=$L/out_units.txt; UL=$L/out_leaves.txt; score $U "level $lv (tree cut $g)"
done
cp -f $U $O/final_units.txt; cp -f $UL $O/final_leaves.txt; log "DONE: $(awk '$1 > 0' $U | sort -n | uniq | wc -l) states"
