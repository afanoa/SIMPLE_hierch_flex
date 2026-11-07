#!/bin/bash
#SBATCH -J klem_merge13
#SBATCH -p norm
#SBATCH -c 32
#SBATCH --mem=240G
#SBATCH -t 16:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# merge-up of an rflex run's leaves on the REAL half-stacks with splitbench13 (junk screen, cross-half EM <= 20 A, re-test beyond the EM
# band with every group map docked on the tested band, stratified null), junk-aware (bu_merge13.py): leaves -> integer ids -> leaf tree
# (state_merge --list) -> per-row poses from each leaf's own node refinement -> per-half pose tables (meta_rows) -> nested tree cuts
# G = N/2 .. 1 -> final units -> unit maps (particle-weighted leaf averages, lp 10, box 160) in view_rflex/<RT>_merged/ + closed-arm readout.
# env RT, PAR (4)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_${RT:?}; P=$R/merge13; mkdir -p $P && cd $P
H=$T/simrnasep692; MRA=${H}_A/meta_rows.txt; MRB=${H}_B/meta_rows.txt; CP=$T/benchcache_rnasep692; log(){ echo "[$(date '+%F %T') merge13 $RT] $*"; }
awk '{print NR, $3, $4, $5}' $R/leaves.txt > $P/leaves_int.txt; NL=$(grep -c . $P/leaves_int.txt); log "$NL leaves"
[ -s $P/tree/merge_tree.txt ] || $PY $T/state_merge.py $P/tree --list $P/leaves_int.txt --link complete --ncut 2 > $P/tree.log 2>&1
[ -s $P/tree/merge_tree.txt ] || { log "leaf tree FAILED"; exit 1; }
if [ ! -s $P/poses_B.txt ]; then
  for pj in $(awk '{print $2}' $R/leaves.txt | xargs -n1 dirname | sort -u); do f=$(ls $pj/*.simple | head -1); k=$(echo $pj | md5sum | cut -c1-12); [ -s $P/dump_$k.txt ] || simple_exec prg=print_project_field projfile=$f oritype=ptcl3D mkdir=no > $P/dump_$k.txt 2> /dev/null; done
  $PY - <<PYEOF
import numpy as np, os, hashlib
L = [l.split() for l in open('$R/leaves.txt') if l.strip()]; out = None
for l in L:
    k = hashlib.md5((os.path.dirname(l[1]) + '\n').encode()).hexdigest()[:12]; d = open('$P/dump_%s.txt' % k).read().splitlines()
    if out is None: out = list(d)
    lab = np.loadtxt(l[5], dtype=int)
    for r in np.where(lab == int(l[6]))[0]: out[r] = d[r]
open('$P/poses_combined.txt', 'w').write('\n'.join(out) + '\n'); print('pose table: %d rows' % len(out))
PYEOF
  for h in A B; do mr=$MRA; [ $h = B ] && mr=$MRB; awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $P/poses_combined.txt $mr > $P/poses_$h.txt; done; fi
log "per-half poses: $(grep -c . $P/poses_A.txt) / $(grep -c . $P/poses_B.txt) rows"
one(){ d=$1; [ -s $d/children_B.txt ] && return 0
  $PY -u $T/splitbench13.py --sim ${H}_A,${H}_B --poses $P/poses_A.txt,$P/poses_B.txt --cand $d/cand_A.txt,$d/cand_B.txt --out $d --cache $CP --ctfabs 1 --pair identity --softem 0 --arch 0 --strata 1 \
    --xem 3 --retest_beyond_em 1 --em_lim 20 --em_min 1000 --dock 2 --junk_proj 0.75 --junk_ssnr 0.35 --summary $d/summary.txt > $d/bench.log 2>&1 \
    && echo "   $(basename $d): $(($(grep -c . $d/labels.txt) - 1)) units -> $(head -1 $d/summary.txt | cut -c1-120)" || echo "   $(basename $d) FAILED: $(grep -a -m1 -A2 'Error\|Traceback' $d/bench.log | tail -1)"; }
U=$R/final_labels.txt; UL=-; JK=-; g=$NL; lv=0
while [ $g -gt 1 ]; do g=$(( (g + 1) / 2 )); lv=$((lv + 1)); L=$P/level$lv
  $PY $T/bu_merge13.py prep $U $UL $JK $P/tree/merge_tree.txt $L $g $MRA $MRB | sed "s/^/level $lv: /"
  while read d us; do case "$d" in \#*) continue;; esac; [ -s $d/cand_A.txt ] || continue
    while [ $(jobs -rp | wc -l) -ge ${PAR:-4} ]; do sleep 10; done; one $d &
  done < $L/nodes.txt; wait
  $PY $T/bu_merge13.py collect $U $UL $JK $L $L/out $MRA $MRB | sed "s/^/level $lv: /"; U=$L/out_units.txt; UL=$L/out_leaves.txt; JK=$L/out_junk.txt
done
cp -f $U $P/final_units.txt; cp -f $UL $P/final_leaves.txt; cp -f $JK $P/final_junk.txt 2>/dev/null; log "$(awk '$1 > 0' $U | sort -n | uniq | wc -l) final units ($(grep -c . $P/final_junk.txt 2>/dev/null || echo 0) junk)"
V=$K_/view_rflex/${RT}_merged; $PY $T/unit_maps.py $V $P/leaves_int.txt $P/final_leaves.txt --junk $P/final_junk.txt
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $P/unit_arm.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | tail -1
awk 'NR == FNR {n[$1]++; next} !/^#/ {k = $1; sub(/.*unit_/, "", k); sub(/\.mrc/, "", k); printf "   %s (%d ptcls): closed-arm %.2f open-arm %.2f %s\n", $1, n[k + 0], $2, $3, $4}' $P/final_units.txt $P/unit_arm.txt | sed 's#[^ ]*/##' | sort -k5,5gr | head -20
log DONE
